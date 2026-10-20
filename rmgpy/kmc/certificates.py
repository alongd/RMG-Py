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
entry could match a long-chain analogue (see its docstring); with a growth
provider the analogues are built at the declared growth sites and compared by
isomorphism.  Anything but an ``excluded:`` verdict (unknown, matched) is
refused (``exclusion_refused``).  The tree estimate (rate rules, averages, auto-generated walk) is
recorded as in ``channels.selection_certificate``.

Spin and domain: every species of the channel carries its spin verdict
(``domain.spin_verdict``) and its v1 domain verdict (``domain.check``); the
multiplicity convention of the bridge is a global guard.

``recertify`` regenerates the mapped reaction on an extended witness, re-runs
the selection and compares through the atom map (``selection_view``).
"""

import hashlib
import json
import re
from collections import OrderedDict, defaultdict
from typing import Dict, Iterable, Optional, Sequence, Tuple

from rmgpy.kmc import channels as chm
from rmgpy.kmc import domain as dom
from rmgpy.kmc import rmg_moves as rm

CERTIFICATE_VERSION = "kmc-rate-certificate/4"
RECORD_CERTIFICATE_SCHEMA = "kmc-record-rate-certificate/3"
REFERENCE_SCHEMA = "kmc-rate-certificate-ref/2"
# record statuses that fire in the simulator (refused records never fire)
FIRING_STATUSES = ("enabled", "irreversible")
CERTIFICATION_MODE = "required"
CERTIFICATION_COVERAGE = "every record"
RADIUS_UNIT = ("ordinary bonds of the witness graph (bond order ignored); the radius of a "
               "group is the largest bond distance from its labelled atoms to any group atom")
H_HANDLING = ("explicit: hydrogens are graph atoms, a bond to H is one bond, an H inside a "
              "ball is listed with its own row")
REPEAT_UNIT_FORMULA = (("C", 8), ("H", 8))
# the domain policy of a certificate without "domain_policy": the module
# domain (``domain.check``), the policy ``windows.LEGACY_V1`` names
LEGACY_DOMAIN_POLICY = "legacy_v1"
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
                if repr(comp) not in self.nodes:
                    self.nodes[repr(comp)] = {
                        "kind": "logic", "operator": type(comp).__name__,
                        "invert": bool(comp.invert),
                        "components": [c if isinstance(c, str) else repr(c)
                                       for c in comp.components]}
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
        by_group = defaultdict(set)
        for m in maps:
            for s, g in m.items():
                by_group[index[g]].add(s.id)
        # ``least`` depends on the id order (audit only); ``by_group`` (every
        # group atom's images over all embeddings) and ``images`` do not and
        # are what support comparison reads
        return {"count": len(keyed), "least": keyed[0],
                "by_group": [[g, sorted(ids)] for g, ids in sorted(by_group.items())],
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


_SELECTION_TABLES = {}


def depository_selection(family, reactants, products) -> Optional[dict]:
    """RMG's depository selection on one reaction (``depository_search``'s
    ``selected``: the entry RMG would rate it from, None for none), searched
    only over the entries whose species formula multisets equal the
    reaction's on one orientation (necessary for ``Reaction.is_isomorphic``),
    so a reaction no entry can match costs no isomorphism test.  The tables
    are kept per family object."""
    from types import SimpleNamespace
    held = _SELECTION_TABLES.get(id(family))
    if held is None or held[0] is not family:
        held = (family, _depository_formulas(family), max_depository_heavy_atoms(family))
        _SELECTION_TABLES[id(family)] = held
    _family, rows, hmax = held
    if any(_heavy(m) > hmax for m in list(reactants) + list(products)):
        return None
    r = _side_formulas(reactants, [{}] * len(reactants))
    p = _side_formulas(products, [{}] * len(products))
    wanted = {(label, index) for label, index, er, ep in rows if (er, ep) in ((r, p), (p, r))}
    if not wanted:
        return None
    stub = SimpleNamespace(depositories=[
        SimpleNamespace(label=d.label, entries={k: e for k, e in d.entries.items()
                                                if (d.label, int(e.index)) in wanted})
        for d in family.depositories])
    return depository_search(stub, reactants, products)["selected"]


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


EXCEPT_MATCHED = "excluded:except_matched"


def exclusion_refused(exclusion) -> bool:
    """Whether a depository exclusion refuses its channel: every verdict but
    an ``excluded:`` argument (``refused:unknown``, none).  An
    ``excluded:except_matched`` channel certifies: its finite matched sizes
    are depository overrides (``exclusion_matches``)."""
    return not (exclusion or {}).get("verdict", "").startswith("excluded:")


def exclusion_matches(exclusion) -> list:
    """The matched analogues an ``excluded:except_matched`` exclusion
    certifies with (each one a depository override), else []."""
    exclusion = exclusion or {}
    return list(exclusion.get("matches") or ()) if exclusion.get("verdict") == EXCEPT_MATCHED \
        else []


def reaction_growth_sites(growth, reactants, products) -> Optional[tuple]:
    """((reactant index, Port, coordinate), ...): the provider's growth sites
    on each reactant (``growth.growth_sites``) that the reaction leaves
    alone (port atom and caps have the same row, bonds included, on both
    sides); None when the provider cannot say for some reactant."""
    before = {a: row for m in reactants for a, row in rm.molecule_rows(m).items()}
    after = {a: row for m in products for a, row in rm.molecule_rows(m).items()}
    out = []
    for k, mol in enumerate(reactants):
        sites = growth.growth_sites(mol)
        if sites is None:
            return None
        for port, coordinate in sites:
            if all(before[a] == after.get(a) for a in (port.atom_id,) + port.cap_atom_ids):
                out.append((k, port, coordinate))
    return tuple(out)


def _insert_unit(side, port, fragment, ids) -> None:
    """Replace ``port``'s caps by ``fragment`` (local index k -> ``ids[k]``)
    in ``side`` (atom id -> [attributes, {neighbour: order}])."""
    for cap in port.cap_atom_ids:
        for other in side.pop(cap)[1]:
            if other in side:
                side[other][1].pop(cap, None)
    for k, attrs in enumerate(fragment.atoms):
        side[ids[k]] = [tuple(attrs), {}]
    bonds = [(ids[i], ids[j], o) for i, j, o in fragment.bonds]
    bonds += [(port.atom_id, ids[i], o) for i, o in fragment.attachments]
    for a, b, o in bonds:
        side[a][1][b] = float(o)
        side[b][1][a] = float(o)


def grown_analogue(growth, reactants, products, sites, vector) -> tuple:
    """(reactants, products) of the analogue that inserts ``vector[i]``
    periods at site ``i`` (``reaction_growth_sites``) on both sides alike,
    each through the provider's ``fragment``; fresh ids above every witness
    id."""
    from dataclasses import replace
    sides = [{a: [row[:5], dict(row[5])] for m in mols for a, row in rm.molecule_rows(m).items()}
             for mols in (reactants, products)]
    fresh = max(max(side) for side in sides) + 1
    for (_k, port, coordinate), count in zip(sites, vector):
        for _ in range(count):
            fragment = growth.fragment(port, coordinate)
            ids = list(range(fresh, fresh + len(fragment.atoms)))
            fresh += len(fragment.atoms)
            for side in sides:
                _insert_unit(side, port, fragment, ids)
            coordinate = next((u, p) for i, u, p in fragment.coordinates
                              if i == fragment.new_port)
            port = replace(port, atom_id=ids[fragment.new_port],
                           cap_atom_ids=tuple(ids[c] for c in fragment.new_caps))
    return tuple(tuple(rm.molecule_from_rows(part) for part in rm.split_rows(
        {a: tuple(attrs) + (tuple(sorted(bonds.items())),) for a, (attrs, bonds) in side.items()}))
        for side in sides)


def _heavy(mol) -> int:
    return sum(1 for a in mol.atoms if not a.is_hydrogen())


def _growth_vectors(sites, reactants, products, unit_heavy, hmax):
    """Every vector of periods per site, at least one in total, that keeps
    every grown species within ``hmax`` heavy atoms."""
    product_of = {a.id: j for j, m in enumerate(products) for a in m.atoms}
    owners = [(k, product_of[port.atom_id]) for k, port, _c in sites]
    sizes = ([_heavy(m) for m in reactants], [_heavy(m) for m in products])
    vector = [0] * len(sites)

    def fits():
        grown = ([0] * len(reactants), [0] * len(products))
        for (k, j), count in zip(owners, vector):
            grown[0][k] += count * unit_heavy
            grown[1][j] += count * unit_heavy
        return all(base + extra <= hmax for side in (0, 1)
                   for base, extra in zip(sizes[side], grown[side]))

    def walk(i):
        if i == len(sites):
            if any(vector):
                yield tuple(vector)
            return
        vector[i] = 0
        while fits():
            yield from walk(i + 1)
            vector[i] += 1
        vector[i] = 0

    yield from walk(0)


def _depository_entries(family) -> dict:
    return {(depository.label, int(entry.index)): entry
            for depository in family.depositories for entry in depository.entries.values()}


def _side_formulas(mols, extra) -> tuple:
    out = []
    for mol, count in zip(mols, extra):
        f = _formula(mol)
        for element, n in count.items():
            f[element] = f.get(element, 0) + n
        out.append(tuple(sorted(f.items())))
    return tuple(sorted(out))


def decide_by_construction(family, reactants, products, growth, unit, hmax, candidates, cache):
    """The analogues of this reaction at the growth provider's sites, each
    compared with every formula candidate by ``Reaction.is_isomorphic`` in
    the orientation(s) its formulas allow (forward: the entry's reactants are
    the analogue's reactants; reverse: its products).  Returns (sites,
    analogues compared, matches); sites None when the provider cannot say."""
    from rmgpy.reaction import Reaction
    sites = reaction_growth_sites(growth, reactants, products)
    if sites is None:
        return None, 0, []
    ekey = ("entries", family.label)
    if ekey not in cache:
        cache[ekey] = _depository_entries(family)
    entries = cache[ekey]
    fkey = ("formulas", family.label)
    rows = {(label, index): (er, ep) for label, index, er, ep in cache[fkey]}
    wanted = [(c["depository"], c["index"]) for c in candidates]
    unit_heavy = sum(v for k, v in unit.items() if k != "H")
    product_of = {a.id: j for j, m in enumerate(products) for a in m.atoms}
    compared, matches, seen = 0, [], set()
    for vector in _growth_vectors(sites, reactants, products, unit_heavy, hmax):
        extra = ([defaultdict(int) for _ in reactants], [defaultdict(int) for _ in products])
        for (k, port, _c), count in zip(sites, vector):
            for element, n in unit.items():
                extra[0][k][element] += count * n
                extra[1][product_of[port.atom_id]][element] += count * n
        grown = (_side_formulas(reactants, extra[0]), _side_formulas(products, extra[1]))
        tests = [(key, orientation) for key in wanted
                 for orientation, sides in (("forward", rows[key]), ("reverse", rows[key][::-1]))
                 if tuple(sides) == grown]
        if not tests:
            continue
        new_r, new_p = grown_analogue(growth, reactants, products, sites, vector)
        forward = Reaction(reactants=[_species(m) for m in new_r],
                           products=[_species(m) for m in new_p])
        reverse = Reaction(reactants=forward.products, products=forward.reactants)
        compared += 1
        for (label, index), orientation in tests:
            rxn = forward if orientation == "forward" else reverse
            if not entries[(label, index)].item.is_isomorphic(rxn, either_direction=False):
                continue
            smiles = (sorted(m.to_smiles() for m in new_r), sorted(m.to_smiles() for m in new_p))
            key = (label, index, orientation, tuple(smiles[0]), tuple(smiles[1]))
            if key in seen:
                continue
            seen.add(key)
            entry = entries[(label, index)]
            matches.append({"depository": label, "index": index,
                            "label": entry.label, "orientation": orientation,
                            "periods": list(vector), "reactants": smiles[0],
                            "products": smiles[1],
                            "species_heavy_atoms": {"reactants": sorted(map(_heavy, new_r)),
                                                    "products": sorted(map(_heavy, new_p))},
                            "rank": entry.rank, "kinetics": repr(getattr(entry, "data", None)),
                            # the analogue itself, family template direction
                            "reaction": {"reactants": sorted(m.to_adjacency_list() for m in new_r),
                                         "products": sorted(m.to_adjacency_list() for m in new_p)}})
    return sites, compared, matches


def depository_exclusion(family, reactants, products, *, cache=None,
                         repeat=REPEAT_UNIT_FORMULA, growth=None) -> dict:
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
    3. Without a ``growth`` provider: ``refused:unknown`` with the candidate
       entries.  With one (the declared chain, ``emit.certification_rules``)
       the analogues are built, not counted: the provider names the growth
       sites of each reactant (``growth_sites``; none at an atom the reaction
       edits), every vector of whole periods at those sites within ``H_max``
       is inserted on both sides alike through the provider's ``fragment``,
       and each analogue whose formulas are a candidate's is compared with it
       by ``Reaction.is_isomorphic`` (``decide_by_construction``).
       ``excluded:isomorphism`` when none matches; ``excluded:except_matched``
       when some do (``matched`` names the entries; ``matches`` lists each
       analogue: periods per site, species and heavy-atom sizes, orientation,
       the entry's rank and kinetics as RMG reads them, the analogue's
       adjacency lists): certified, RMG rates those finite sizes from the
       depository and the catalog overrides them there; ``refused:unknown``
       when the provider names no growth sites or grows by another unit."""
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
    if growth is None:
        out.update({"verdict": "refused:unknown", "candidates": candidates,
                    "argument": "a depository entry has the formulas of a grown witness; "
                                "whether it is isomorphic to an actual long-chain analogue is "
                                "not decided"})
        return out
    out.update({"candidates": candidates, "growth_provider": growth.name})
    if dict(growth.unit_formula() or ()) != rep:
        out.update({"verdict": "refused:unknown",
                    "argument": f"the growth provider {growth.name} grows by "
                                f"{dict(growth.unit_formula() or ())}, not by the witness's "
                                f"growth unit {rep}: the analogues are not decided"})
        return out
    sites, compared, matches = decide_by_construction(
        family, reactants, products, growth, rep, hmax, candidates, cache)
    if sites is None:
        out.update({"verdict": "refused:unknown",
                    "argument": f"the growth provider {growth.name} names no growth sites for "
                                "a witness species: where a unit inserts is not known, so the "
                                "analogues are not decided"})
        return out
    out.update({"growth_sites": [{"reactant": k, "atom": port.atom_id, "tag": port.tag}
                                 for k, port, _c in sites],
                "analogues_compared": compared, "matches": matches})
    if not matches:
        out.update({"verdict": "excluded:isomorphism",
                    "argument": "every long-chain analogue (whole periods inserted at the "
                                "declared growth sites, on both sides, within the largest "
                                "depository species) whose formulas are a candidate's was "
                                "compared with it by Reaction.is_isomorphic in both "
                                "orientations: none is isomorphic"})
        return out
    named = sorted({(m["depository"], m["index"]) for m in matches})
    out.update({"verdict": EXCEPT_MATCHED,
                "matched": [f"{d}#{i}" for d, i in named],
                "argument": "every long-chain analogue was compared as for excluded:isomorphism; "
                            "the analogues in matches (each a finite size: species within the "
                            "largest depository species) ARE depository entries, so RMG rates "
                            "exactly those oligomers from the depository (entry kinetics as RMG "
                            "reads them in matches) and every other analogue by the rule "
                            "estimate: certified, with the catalog carrying each match as a "
                            "size-conditioned override (template_emission.depository_overrides)"})
    return out


# ---------------------------------------------------------------------------
# bounded caches
# ---------------------------------------------------------------------------
class _Bounded:
    """A least-recently-used mapping with a hard size bound."""

    def __init__(self, maxsize: int):
        if int(maxsize) < 1:
            raise ValueError("a cache bound must be positive")
        self.maxsize = int(maxsize)
        self._data = OrderedDict()

    def __contains__(self, key):
        return key in self._data

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        self._data.move_to_end(key)
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value
        self._data.move_to_end(key)
        while len(self._data) > self.maxsize:
            self._data.popitem(last=False)

    def get(self, key, default=None):
        if key in self._data:
            return self[key]
        return default

    def values(self):
        return list(self._data.values())


# ---------------------------------------------------------------------------
# mapped moves: identity, edits, matched form
# ---------------------------------------------------------------------------
def _row_edges(row):
    return {int(n): float(o) for n, o in row[5]}


def move_edits(reactant_rows, product_rows) -> dict:
    """The recipe's edits of one move: changed atom attributes (element,
    isotope, charge, radical electrons, lone pairs) and changed bonds,
    between the reactant forms the move used and its product (one id
    namespace; explicit H)."""
    atoms, bonds = [], []
    for a in sorted(reactant_rows):
        before, after = reactant_rows[a], product_rows.get(a)
        if after is None:
            raise CertificateError(f"atom {a} is missing from the move's product")
        if tuple(before[:5]) != tuple(after[:5]):
            atoms.append([a, list(before[:5]), list(after[:5])])
    pairs = set()
    for rows in (reactant_rows, product_rows):
        for a, row in rows.items():
            for n, _o in row[5]:
                pairs.add((min(a, int(n)), max(a, int(n))))
    for a, b in sorted(pairs):
        o1 = _row_edges(reactant_rows[a]).get(b, 0.0)
        o2 = _row_edges(product_rows[a]).get(b, 0.0)
        if o1 != o2:
            bonds.append([a, b, o1, o2])
    return {"atoms": atoms, "bonds": bonds}


def _form_bonds(reactant_rows, mapped_ids) -> list:
    """Bonds (with their orders in the matched resonance forms) among the
    atoms the template subgraph isomorphism mapped."""
    out = []
    for a in sorted(mapped_ids):
        for n, o in reactant_rows[a][5]:
            if n in mapped_ids and a < n:
                out.append([a, int(n), float(o)])
    return out


def mapping_key(mapping: dict, t=None) -> tuple:
    """Identity of one certified mapping, ids through ``t``: the template
    slot of every labelled atom, the full ``Move.mapped`` correspondence
    (slot, template structure, template atom, atom), the recipe's edits and
    the bond orders of the matched resonance forms among the mapped atoms."""
    t = t or (lambda a: a)
    labels = tuple(sorted((int(s), lab, t(a)) for s, lab, a in mapping["labels"]))
    mapped = tuple(sorted((int(s), int(i), int(j), t(a))
                          for s, i, j, a in mapping.get("mapped", [])))
    edits = mapping.get("edits") or {"atoms": [], "bonds": []}
    e_atoms = tuple(sorted((t(a), tuple(b), tuple(c)) for a, b, c in edits["atoms"]))
    e_bonds = tuple(sorted((tuple(sorted((t(a), t(b)), key=repr)), o1, o2)
                           for a, b, o1, o2 in edits["bonds"]))
    forms = tuple(sorted((tuple(sorted((t(a), t(b)), key=repr)), o)
                         for a, b, o in mapping.get("form_bonds", [])))
    return labels, mapped, e_atoms, e_bonds, forms


def _transport_fn(transport):
    if transport is None:
        return lambda a: a
    return lambda a: transport.get(a, ("unmapped", a))


# ---------------------------------------------------------------------------
# certifier
# ---------------------------------------------------------------------------
class Certifier:
    """Builds certificates.  Every cache is bounded (least recently used):
    channel sets (``max_channel_sets``; the resonance forms of a channel set
    live and leave with it), descents with their node semantics, depository
    searches, record certificates, RMG selection checks and per-family
    tables (``max_entries`` each).

    ``domain`` and ``repeat``: the species domain and the depository growth
    unit of RECORD certification (``certify_record``) and of record-level
    recertification (``recertify``) when none is given there; None keeps the
    module domain and ``REPEAT_UNIT_FORMULA`` (an artifact compiled from a
    declaration passes the declaration's own: ``emit.certification_rules``).
    ``growth``: the provider whose growth sites decide every depository
    exclusion by construction (``depository_exclusion``; its unit is the
    growth unit where none is given); None keeps the formula argument."""

    def __init__(self, families, resonance: str = "all", max_channel_sets: int = 64,
                 max_entries: int = 4096, domain=None, repeat=None, growth=None):
        self.families = families
        self.domain = domain
        self.repeat = repeat
        self.growth = growth
        self.resonance = resonance
        self.max_channel_sets = max_channel_sets
        self._channel_sets = _Bounded(max_channel_sets)
        self._form_sets = _Bounded(max_channel_sets)
        self._descents = _Bounded(max_entries)
        self._depository = _Bounded(max_entries)
        self._records = _Bounded(max_entries)
        self._rmg = _Bounded(max_entries)
        self.cache = _Bounded(max_entries)     # per-family depository tables

    def cache_sizes(self) -> dict:
        return {"channel_sets": len(self._channel_sets), "forms": len(self._form_sets),
                "descents": len(self._descents), "depository": len(self._depository),
                "records": len(self._records), "rmg_checks": len(self._rmg),
                "family_tables": len(self.cache)}

    # -- forms and labelled structures --------------------------------------
    def channel_set(self, family, reactants):
        """Forward channels on the aromatic roots of ``reactants`` (cached by
        the witness's integer-keyed signature; bounded)."""
        roots = [aromatic_root(m) for m in reactants]
        key = ("channels", family.label, tuple(rm._signature(rm.molecule_rows(m))
                                                for m in roots))
        cs = self._channel_sets.get(key)
        if cs is None:
            cs = chm.build_channels(family, roots, "forward", resonance=self.resonance)
            self._channel_sets[key] = cs
        return cs

    def _forms(self, channel_set, k):
        entry = self._form_sets.get(id(channel_set))
        if entry is None or entry[0] is not channel_set:
            entry = (channel_set, {})
            self._form_sets[id(channel_set)] = entry
        forms = entry[1]
        if k not in forms:
            forms[k] = rm.reactant_forms(channel_set.reactants[k], channel_set.resonance)
        return forms[k]

    def _reactant_rows(self, channel_set, move):
        rows = {}
        for k in range(len(channel_set.reactants)):
            rows.update(rm.molecule_rows(self._forms(channel_set, k)[move.forms[k]]))
        return rows

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
        # template order is kept: the descents record participant indices
        return (family.label, tuple(
            (rm._signature(rm.molecule_rows(m)),
             tuple(sorted((a.id, a.label) for a in m.atoms if a.label))) for m in molecules))

    def descents(self, family, channel_set, move, nodes, store):
        """Ref of the recorded descents of one mapping (template direction);
        the semantics of every node they reference are copied into
        ``nodes`` (also on a cache hit)."""
        mols = self._labelled(channel_set, move)
        key = self._descent_key(family, mols)
        cached = self._descents.get(key)
        if cached is None:
            local = {}
            cached = (template_descents(family, mols, local), local)
            self._descents[key] = cached
        descents, semantics = cached
        nodes.update(semantics)
        ref = _digest(key)
        store[ref] = descents
        return ref

    # -- channel ------------------------------------------------------------
    def certify(self, family, reactants, *, pathway=None, products=None, member_labels=None,
                channel_set=None, species=None, repeat=None, domain=None) -> dict:
        """Certificate of the forward (family-template direction) channel(s) on
        ``reactants`` selected by ``pathway`` + product ``products`` (species
        keys), or by ``member_labels`` (a set of frozensets of (label, atom id):
        every channel with a member one of whose mappings has those labels).
        ``repeat``: the growth unit of ``depository_exclusion`` (element,
        count) pairs; None keeps ``REPEAT_UNIT_FORMULA``.  ``domain``: the
        species domain of the channel eligibility (an object with ``name``,
        ``check`` and ``spin_verdict``, ``windows.PolicyDomain``); None is the
        module domain (``domain.check``) and leaves the certificate as it
        always was, a domain adds ``"domain_policy": domain.name``."""
        species = species or chm.SpeciesKeys()
        given = list(reactants) if channel_set is None else list(channel_set.reactants)
        input_refusals = []
        for k, m in enumerate(given):
            unpaired = sum(int(a.radical_electrons) for a in m.atoms)
            if int(m.multiplicity) != unpaired + 1:
                input_refusals.append(
                    f"input reactant {k}: multiplicity {m.multiplicity} with {unpaired} "
                    "unpaired electrons violates the declared maximum-multiplicity convention")
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
                                          for m in channel_set.reactants],
                            "atoms": _witness_atoms(channel_set.reactants)},
                "channels": [], "descents": {}, "nodes": {}, "species": {}, "gaps": [],
                "input_refusals": input_refusals}
        if not chosen:
            cert["gaps"].append("no channel on the witness matches the selector")
        for ch in chosen:
            cert["channels"].append(self._channel(family, channel_set, ch, cert, species,
                                                  repeat=repeat, domain=domain))
        cert["global_guards"] = _global_guards(cert)
        cert["gaps"].extend(completeness_gaps(cert))
        cert["complete"] = not cert["gaps"]
        # completeness is about what RMG read; whether a long-chain analogue
        # could switch the rate to a depository entry is a separate verdict
        cert["depository_exclusion_refused"] = [
            ";".join(ch["pathway"]) for ch in cert["channels"]
            if not ch.get("whole_species") and exclusion_refused(ch.get("depository_exclusion"))]
        # domain eligibility: its own disposition, consumed by record
        # certification (certified) and by recertification
        cert["eligible"] = bool(cert["channels"]) and not input_refusals and all(
            (ch.get("eligibility") or {}).get("disposition") == "eligible"
            for ch in cert["channels"])
        if domain is not None:
            cert["domain_policy"] = domain.name
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

    def _moves_by_outcome(self, channel_set, ch):
        by_outcome = defaultdict(list)
        for m in self._member_moves(channel_set, ch):
            by_outcome[rm.product_key(m)].append(m)
        return {o: sorted(v, key=rm.Move.sort_key) for o, v in by_outcome.items()}

    def _channel(self, family, channel_set, ch, cert, species, repeat=None, domain=None):
        out = {"direction": ch.direction, "pathway": list(ch.pathway), "kind": ch.kind,
               "reactants": list(ch.reactants), "products": list(ch.products),
               "size": ch.size, "mappings": list(ch.mappings), "support": ch.support,
               "ambiguous_with": [list(p) for p in ch.ambiguous_with]}
        if ch.direction != "forward":
            out["gap"] = "certificates are built in the family template (forward) direction"
            out["eligibility"] = {"disposition": "refused:direction", "refused_species": {}}
            out["publishable"] = False
            return out
        by_outcome = self._moves_by_outcome(channel_set, ch)
        others = channel_set.by_outcome()
        members = []
        for outcome in ch.members:
            moves = by_outcome.get(outcome, [])
            maps = []
            for m in moves:
                rows = self._reactant_rows(channel_set, m)
                mapped_ids = {int(a) for _s, _i, _j, a in m.mapped}
                identity = {"order": list(m.order), "forms": list(m.forms),
                            "labels": [list(x) for x in m.labels],
                            "mapped": [list(x) for x in m.mapped],
                            "edits": move_edits(rows, m.product_rows()),
                            "form_bonds": _form_bonds(rows, mapped_ids), "status": m.status}
                try:
                    ref = self.descents(family, channel_set, m, cert["nodes"], cert["descents"])
                except chm.UnsupportedSelection as exc:
                    maps.append({**identity, "descent": None, "gap": str(exc)})
                    continue
                maps.append({**identity, "descent": ref,
                             "leaves": descent_leaves(cert["descents"][ref])})
            amb = sorted({c.pathway for c in others.get(outcome, ()) if c.pathway != ch.pathway})
            members.append({"outcome": _digest(outcome), "mappings": maps,
                            "ambiguous_with": [list(p) for p in amb]})
        out["members"] = members
        # species verdicts and the channel's domain eligibility
        rep = ch.representatives[0]
        mols = list(channel_set.reactants) + [rm.molecule_from_rows(p) for p in
                                              rm.split_rows(rep.product_rows())]
        refused = {}
        for mol in mols:
            k = species(mol)
            if k not in cert["species"]:
                if domain is None:
                    spin, verdict = dom.spin_verdict(mol), dom.check(mol)
                else:
                    spin, verdict = domain.spin_verdict(mol), domain.check(mol)
                cert["species"][k] = {
                    "spin": {"status": spin.status, "reason": spin.reason,
                             "radical_electrons": spin.radical_electrons,
                             "convention": spin.convention},
                    "domain": verdict.as_dict()}
            if not cert["species"][k]["domain"]["accepted"]:
                refused[k] = cert["species"][k]["domain"]["reason"]
        out["eligibility"] = {"disposition": "refused:domain" if refused else "eligible",
                              "refused_species": refused}
        out["publishable"] = ch.support == chm.SUPPORTED and not refused
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
            if repeat is None:
                repeat = REPEAT_UNIT_FORMULA if self.growth is None else \
                    self.growth.unit_formula()
            out["depository_exclusion"] = depository_exclusion(
                family, reactants, products, cache=self.cache, repeat=repeat,
                growth=self.growth)
            try:
                kind, sel = self._tree_selection(family, ch)
                out["source_kind"], out["selection"] = kind, sel
            except chm.UnsupportedSelection as exc:
                out["source_kind"], out["selection"] = "unsupported", {"reason": str(exc)}
        out["rmg_check"] = self._rmg_check(family, channel_set, ch, out)
        return out

    def _tree_selection(self, family, ch):
        """The tree estimate RMG makes for the channel's template.  For an
        auto-generated tree the uncertainty walk is evaluated first and a
        walk RMG would never leave is refused BEFORE get_kinetics_for_template
        is called (RMG's own loop spins forever there)."""
        template = family.retrieve_template(list(ch.pathway))
        if family.auto_generated:
            walk = chm.autogenerated_walk(family, ch.pathway)
            if not walk["rmg_terminates"]:
                chm._autogenerated_selection(family, ch.pathway, None, walk=walk)  # raises
                raise chm.UnsupportedSelection("RMG's estimate_kinetics walk does not terminate")
            kinetics, entry = family.get_kinetics_for_template(template, degeneracy=1)
            sel = chm._autogenerated_selection(family, ch.pathway, kinetics, walk=walk)
            sel["rmg_entry"] = str(entry) if entry else None
            return "autogenerated_tree", sel
        kinetics, entry = family.get_kinetics_for_template(template, degeneracy=1)
        kind, sel = chm._rate_rule_selection(family, ch.pathway, kinetics)
        sel["supplying_node"] = sel["supplying_templates"][0] \
            if len(sel["supplying_templates"]) == 1 else None
        sel["rmg_entry"] = str(entry) if entry else None
        return kind, sel

    def _rmg_selection(self, family, channel_set, ch, move) -> dict:
        """RMG's own production call on one mapped move: ``get_kinetics``
        (template labels, degeneracy 1, ``return_all_kinetics=False``, as the
        compiler calls it), the selection direction (``is_forward``), the
        rate's uncertainty, for an auto-generated tree the SUPPLYING node
        (RMG's "Estimated from node" comment; the returned entry is None for
        a non-exact selection) and, for a hand-built tree, the compiler's own
        source extraction (``compiler._rate_rule_source``: rules, training
        entries, weights) on the reaction normalized exactly as production
        normalizes it: swapped when RMG selected the rate for the opposite
        direction (``is_forward`` False).  Only called when the selection is
        supported (a non-terminating auto-generated walk is refused
        before)."""
        key = (family.label, tuple(ch.pathway),
               tuple(rm._signature(rm.molecule_rows(m)) for m in channel_set.reactants),
               move.sort_key())
        cached = self._rmg.get(key)
        if cached is not None:
            return cached
        rxn, _r, _p = chm.template_reaction(family, channel_set, move)
        rxn.reactants = [_species(s.molecule[0]) for s in rxn.reactants]
        rxn.products = [_species(s.molecule[0]) for s in rxn.products]
        kinetics, source, entry, is_forward = family.get_kinetics(
            rxn, template_labels=list(ch.pathway), degeneracy=1, return_all_kinetics=False)
        rxn.kinetics = kinetics
        if not is_forward:
            # production (_build_linked_family_pair) swaps the estimate before
            # extracting its contributors: a training entry matched for the
            # opposite direction is then read as reverse=False
            rxn.reactants, rxn.products = rxn.products, rxn.reactants
        tree = isinstance(source, str)
        out = {"source": str(source), "entry": str(entry) if entry else None,
               "rank": getattr(entry, "rank", None),
               "depository_index": None if tree else int(entry.index),
               "is_forward": bool(is_forward), "sources": None,
               "autogenerated": bool(family.auto_generated),
               "uncertainty": str(getattr(kinetics, "uncertainty", None) or "UNKNOWN"),
               "supplier": None}
        if family.auto_generated and tree:
            out["supplier"] = _estimated_from(getattr(kinetics, "comment", "") or "")
        if not family.auto_generated:
            out["sources"] = production_sources(family, rxn)
        self._rmg[key] = out
        return out

    def _rmg_check(self, family, channel_set, ch, out, move=None):
        """RMG's own get_kinetics (on ``move``, default the representative)
        must agree with the recorded selection: the depository entry, and for
        a tree estimate the returned entry and the extracted contributors."""
        if out.get("source_kind") in (None, "unsupported"):
            return {"skipped": "selection unsupported"}
        rmg = self._rmg_selection(family, channel_set, ch, move or ch.representatives[0])
        ours = out["depository"]["selected"]
        problems = []
        if (rmg["depository_index"] is None) != (ours is None) or (
                ours is not None and rmg["depository_index"] != ours["index"]):
            problems.append(f"RMG selected depository entry {rmg['depository_index']}, the "
                            f"recorded search {ours}")
        if ours is None:
            sel = out.get("selection") or {}
            if sel.get("rmg_entry") != rmg["entry"]:
                problems.append(f"RMG returned entry {rmg['entry']!r}, the certified tree "
                                f"selection {sel.get('rmg_entry')!r}")
            if not family.auto_generated and sel.get("sources") != rmg["sources"]:
                problems.append("RMG's contributors differ from the certified contributors")
        if problems and move is None:
            out.setdefault("gaps", []).extend(problems)
        return {"rmg_depository_entry": rmg["depository_index"], "agrees": not problems,
                "problems": problems}

    # -- records ------------------------------------------------------------
    def certify_record(self, record, rate_owner=None) -> dict:
        """Record certificate: the channel certificate of the record's
        template side, bound to the record (``owner``) and to the record's
        production selection (``binding``).

        The rate-owning record is the record itself for an "RMG family
        estimate" and the linked forward record (``rate_owner``, which must
        be its ``reverse_of`` partner) for a "reference-thermo reverse"; a
        library-derived rate (k_library or k_library / Kc) is never an RMG
        selection.  ``binding`` requires (a) members of the certified channel
        whose mapped products equal the record's products under the record's
        heavy-atom map (through the padding projection when padded) and (b)
        RMG's production selection on such a member (source, entry, rank,
        exact flag, rules and training entries with weights, the selection
        direction, the supplying node of an auto-generated tree) equal to the
        rate-owning record's serialized rate source, and (c) the record's
        executable rewrite: a complete canonical heavy-atom map, and its
        complete executable operation stream (``bond_ops`` and
        ``feature_ops``) replayed on its reactant graphs with the runtime's
        semantics reproduces its declared products exactly under that map,
        before any chemistry normalisation (``rewrite_problems``).  A mismatch refuses the binding;
        nothing is re-estimated in its place.  ``certified`` = complete and
        eligible and bound and no refused depository exclusion and every
        bound channel publishable (supported spin/pathway status and domain
        eligible)."""
        record = _as_dict(record)
        rs = record.get("rate_source") or {}
        kind = rs.get("kind")
        owner = record_owner(record)
        head = {"schema": RECORD_CERTIFICATE_SCHEMA, "owner": owner,
                "event_id": record.get("event_id"), "family": record.get("family"),
                "template": record.get("template"), "rate_source_kind": kind,
                "orientation": record.get("orientation")}

        def not_rmg(reason):
            owner.update({"rate_owner_event_id": None, "rate_source_sha256": None})
            return {**head, "certificate": None, "source_kind": "not_rmg_selection",
                    "complete": False, "eligible": False, "certified": False,
                    "binding": {"status": "not_applicable", "reason": reason, "members": [],
                                "mismatches": []},
                    "gaps": [reason]}

        problems = []
        prs = None
        owner_orientation = None
        if kind == "RMG family estimate":
            prs = rs
            owner_orientation = record.get("orientation", "as_generated")
            owner.update({"rate_owner_event_id": record.get("event_id"),
                          "rate_source_sha256": _sha(rs)})
        elif kind == "reference-thermo reverse":
            forward_kind = (rs.get("forward_rate_source") or {}).get("kind")
            if rs.get("reference") == "k_library / Kc" or forward_kind == "kMC kinetics library":
                return not_rmg("library-derived reverse rate (k_library / Kc): not an RMG "
                               "selection")
            owner.update({"rate_owner_event_id": None, "rate_source_sha256": None})
            if rate_owner is None:
                problems.append("rate owner: the linked rate-owning record was not supplied")
            else:
                ro = _as_dict(rate_owner)
                rok = (ro.get("rate_source") or {}).get("kind")
                if rok == "kMC kinetics library":
                    return not_rmg("library-derived reverse rate (its rate-owning record is a "
                                   "kMC kinetics library rate): not an RMG selection")
                if rok != "RMG family estimate":
                    return not_rmg(f"rate-owning record source {rok!r} is not an RMG family "
                                   "selection")
                if ro.get("event_id") != record.get("reverse_of") or \
                        ro.get("reverse_of") != record.get("event_id"):
                    problems.append("rate owner: the supplied record is not this record's "
                                    "reverse_of partner")
                else:
                    prs = ro["rate_source"]
                    owner_orientation = ro.get("orientation", "as_generated")
                    owner.update({"rate_owner_event_id": ro.get("event_id"),
                                  "rate_source_sha256": _sha(prs)})
        else:
            return not_rmg(f"rate source {kind!r} is not an RMG family selection (library or "
                           "archived rate)")
        family = self.families[record["family"]]
        side, mols, others = self.template_side(record)
        pathway = tuple(record["template"].split(";")) if record["template"] else None
        rkey = (record["family"], record["template"],
                tuple(record[side]), tuple(_graph_field(record, side, other=True)))
        cert = self._records.get(rkey)
        if cert is None:
            species = chm.SpeciesKeys()
            products = tuple(sorted(species(m) for m in others))
            cert = self.certify(family, mols, pathway=pathway, products=products,
                                species=species, repeat=self.repeat, domain=self.domain)
            self._records[rkey] = cert
        binding = self._bind(record, family, side, mols, others, cert, prs, problems,
                              owner_orientation)
        padded = (record.get("proxy_padding") or {}).get("status") == "padded"
        gaps = list(cert["gaps"]) + [f"binding: {m}" for m in binding["mismatches"]]
        certified = derived_certified(cert, binding)
        if binding["status"] == "bound" and not certified:
            gaps.extend(f"not publishable: channel {i} ({';'.join(c['pathway'])}): support "
                        f"{c.get('support')}, {(c.get('eligibility') or {}).get('disposition')}"
                        for i, c in _bound_channels(cert, binding) if not c.get("publishable"))
        return {**head, "template_side": side, "padded": padded,
                "role": "rate_source" if kind == "RMG family estimate" else "reverse_of_rate_source",
                "certificate": cert, "complete": cert["complete"], "eligible": cert["eligible"],
                "binding": binding, "certified": certified, "gaps": gaps}

    def _bind(self, record, family, side, mols, others, cert, prs, problems,
              owner_orientation=None) -> dict:
        mismatches = list(problems)
        members = []
        try:
            pins = record_pins(record, side, mols, others)
        except CertificateError as exc:
            pins = None
            mismatches.append(f"member: {exc}")
        if pins is not None:
            mismatches.extend(f"rewrite: {p}" for p in rewrite_problems(record))
        channel_set = self.channel_set(family, mols) if cert["channels"] else None
        matched_move = None
        if pins is not None and channel_set is not None:
            owner_of = {a.id: k for k, m in enumerate(others) for a in m.atoms}
            for ci, cch in enumerate(cert["channels"]):
                ch = _find_channel(channel_set, cch)
                if ch is None:
                    continue
                by_outcome = self._moves_by_outcome(channel_set, ch)
                for k, outcome in enumerate(ch.members):
                    hits = [mi for mi, m in enumerate(by_outcome.get(outcome, []))
                            if _move_matches(m, others, owner_of, pins)]
                    if hits:
                        members.append({"channel": ci, "member": k, "mappings": hits})
                        if matched_move is None:
                            matched_move = (cch, ch, by_outcome[outcome][hits[0]])
            if not members:
                mismatches.append("member: no certified member corresponds to the record's "
                                  "mapped event (its products under the record's atom map)")
        supplier, direction = None, None
        if matched_move is not None and prs is not None:
            cch, ch, move = matched_move
            if cch.get("source_kind") in (None, "unsupported"):
                mismatches.append("selection: the certified selection is unsupported")
            else:
                check = self._rmg_check(family, channel_set, ch, cch, move=move)
                mismatches.extend(f"selection: {p}" for p in check["problems"])
                rmg = self._rmg_selection(family, channel_set, ch, move)
                if cch.get("source_kind") == "autogenerated_tree":
                    supplier = (cch.get("selection") or {}).get("supplying_node")
                direction = "forward" if rmg["is_forward"] else "reverse"
                mismatches.extend(production_mismatches(prs, rmg, ";".join(ch.pathway),
                                                        supplier=supplier))
                # the rate owner (the record itself for a family estimate, its
                # reverse_of partner otherwise) runs in the direction RMG
                # selected the rate for: as_generated iff is_forward
                if (owner_orientation == "as_generated") != bool(rmg["is_forward"]):
                    mismatches.append("direction: rate owner runs opposite to selected "
                                      "kinetics")
        elif prs is None and not problems:
            mismatches.append("rate owner: no production selection to bind")
        if not members and not mismatches:
            mismatches.append("member: nothing was bound")
        return {"status": "bound" if not mismatches else "refused", "members": members,
                "mismatches": mismatches, "supplier": supplier, "direction": direction,
                "evidence": binding_evidence(record),
                "compared": ["source", "entry", "rank", "template", "direction", "supplier",
                             "exact", "rules", "training", "heavy map", "rewrite"]}

    def template_side(self, record):
        """(field name, family-template reactant molecules, other side): the
        as_generated record's reactants, the reversed record's products (the
        rate witness graphs when padded)."""
        record = _as_dict(record)
        a, b = _template_fields(record)
        mols = witness_molecules(record[a])
        others = witness_molecules(record[b], offset=sum(len(m.atoms) for m in mols))
        return a, mols, others

    def recertify(self, record, witness, atom_map, certificate=None, rate_owner=None,
                  domain=None) -> dict:
        """Re-run the selection of a record on an extended witness.

        ``witness``: molecules of the record's family-template reactants on
        the extension (for a reversed record: the extension of its products,
        the side the rate is selected on), with atom ids; ``atom_map``: record
        certificate witness id -> extended witness id.  ``certificate``: the
        record's RECORD certificate (default: certified now, through
        ``rate_owner`` for a reference-thermo reverse).  Refused unless it is
        a record certificate owned by ``record`` and by the supplied
        ``rate_owner`` (the record itself for an RMG family estimate),
        re-derived on this certifier's database, successfully bound and
        certified (``record_certificate_refusals``: the artifact check, one
        code path); its inner certificate is never unwrapped otherwise.  See
        ``recertify_certificate`` (``domain``: the new certificate's domain,
        default the certifier's; the growth unit is the certifier's; a change
        of domain policy is refused)."""
        record = _as_dict(record)
        old = certificate if certificate is not None else \
            self.certify_record(record, rate_owner=rate_owner)
        refusals = record_certificate_refusals(record, old, rate_owner=rate_owner,
                                               certifier=self)
        if refusals:
            return _refused_recertification(refusals)
        return self.recertify_certificate(self.families[record["family"]], old["certificate"],
                                          witness, atom_map, repeat=self.repeat,
                                          domain=self.domain if domain is None else domain)

    def recertify_certificate(self, family, old_cert, witness, atom_map, repeat=None,
                              domain=None) -> dict:
        """Regenerate the mapped reaction on ``witness``, re-run the selection
        and compare through the map.

        Refused (``same_selection`` False, ``refusals`` say why) unless the
        old certificate has at least one certified mapping, is complete,
        domain-eligible and without a refused depository exclusion; the atom
        map transports the whole old witness (injective, elements, bonds and
        radical electrons of every surviving atom; only hydrogens outside
        every certified mapping may be dropped); and the new certificate is
        complete, eligible and without a refused depository exclusion.
        Otherwise ``compare_certificates``: ``same_selection`` (every
        certified mapping, by its full identity, has the same selection and
        descent verdicts), ``same_support`` (the transported support is
        equal) and ``fixed_point`` (both).  ``repeat``: the new certificate's
        depository growth unit (``certify``).  ``domain``: the new
        certificate's species domain (``certify``); refused
        (``refused:domain_policy_changed``) when it is not the old
        certificate's ``domain_policy`` (absent: ``LEGACY_DOMAIN_POLICY``)."""
        transport = {int(k): int(v) for k, v in atom_map.items()}
        before = old_cert.get("domain_policy", LEGACY_DOMAIN_POLICY)
        after = LEGACY_DOMAIN_POLICY if domain is None else domain.name
        if before != after:
            return _refused_recertification(
                [f"refused:domain_policy_changed: {before} -> {after}"])
        refusals = old_certificate_refusals(old_cert)
        if refusals:
            return _refused_recertification(refusals)
        problems = transport_problems(old_cert, witness, transport)
        if problems:
            return _refused_recertification([f"atom map: {p}" for p in problems])
        wanted = set()
        for ch in old_cert["channels"]:
            for mem in ch.get("members", []):
                for mp in mem["mappings"]:
                    wanted.add(frozenset((lab, transport[a]) for _s, lab, a in mp["labels"]))
        new = self.certify(family, witness, member_labels=wanted, repeat=repeat, domain=domain)
        refusals = [f"extension: {r}" for r in old_certificate_refusals(new)]
        comparison = compare_certificates(old_cert, new, transport)
        if refusals:
            return {**comparison, "certificate": new, "same_selection": False,
                    "same_support": False, "same_cohort": False, "fixed_point": False,
                    "refusals": refusals, "reason": "; ".join(refusals)}
        return {**comparison, "certificate": new, "refusals": []}

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


def _witness_atoms(molecules) -> list:
    """[id, element, radical electrons, sorted neighbour ids, isotope, charge,
    lone pairs, sorted [neighbour id, bond order]] per witness atom."""
    out = []
    for m in molecules:
        for a in m.atoms:
            out.append([a.id, a.element.symbol, int(a.radical_electrons),
                        sorted(b.id for b in a.bonds), rm._isotope(a), int(a.charge),
                        int(a.lone_pairs),
                        sorted([b.id, float(bond.order)] for b, bond in a.bonds.items())])
    return sorted(out)


def _constitution(rows) -> list:
    """The resonance-independent part of witness atom rows: id, element,
    isotope and neighbour ids."""
    return sorted([row[0], row[1], row[4] if len(row) > 4 else None, list(row[3])]
                  for row in rows)


def _sha(value) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def production_sources(family, reaction) -> dict:
    """The contributors the compiler serializes for a hand-built-tree rate
    (``compiler._rate_rule_source``: RMG's ``extract_source_from_comments``),
    restricted to the selection: exact flag, rules and training entries with
    their weights."""
    from rmgpy.kmc.compiler import _rate_rule_source
    source = _rate_rule_source(family, reaction)
    return {"exact": source.get("exact"), "rules": source.get("rules"),
            "training": source.get("training")}


def _estimated_from(comment: str):
    """The supplying node of an auto-generated tree estimate, from RMG's
    ``estimate_kinetics`` comment ("Estimated from node X in family F.")."""
    match = re.search(r"Estimated from node (.+?) in family ", comment)
    return match.group(1) if match else None


def stated_supplier(uncertainty):
    """The supplying node a record's serialized rate states: production
    serializes the selected rule's ``RateUncertainty``, whose correlation is
    the rule's node label (auto-generated trees)."""
    if not isinstance(uncertainty, str):
        return None
    match = re.search(r"correlation='([^']*)'", uncertainty)
    return match.group(1) if match else None


def production_mismatches(prs, rmg, template, supplier=None) -> list:
    """Differences between a record's serialized production rate source and
    RMG's selection on the record's mapped member (see certify_record):
    availability; source, entry and rank; template; the selection direction
    (``family_template_direction`` against RMG's ``is_forward``); for an
    auto-generated tree estimate the supplying node (the record's stated
    supplier, RMG's "Estimated from node" and the certified walk's
    ``supplier`` must be one node); for a hand-built tree the contributors
    with their weights."""
    out = []
    if not prs.get("available", False):
        out.append("available: the record's rate is not available")
    for key in ("source", "entry", "rank"):
        if prs.get(key) != rmg.get(key):
            out.append(f"{key}: record {prs.get(key)!r} != RMG {rmg.get(key)!r}")
    if "template" in prs and prs["template"] != template:
        out.append(f"template: record {prs['template']!r} != certified {template!r}")
    direction = "forward" if rmg.get("is_forward", True) else "reverse"
    if prs.get("family_template_direction") is None:
        out.append("direction: the record does not state the direction RMG selected its "
                   f"rate for (RMG: {direction})")
    elif prs["family_template_direction"] != direction:
        out.append(f"direction: record {prs['family_template_direction']!r} != RMG "
                   f"{direction!r}")
    autogenerated_tree = rmg.get("depository_index") is None and (
        rmg.get("autogenerated") or rmg.get("supplier") is not None or supplier is not None)
    if autogenerated_tree:
        stated = stated_supplier(prs.get("uncertainty"))
        if rmg.get("supplier") is None:
            out.append("supplier: RMG's estimate does not name its supplying node")
        if stated is None:
            out.append("supplier: the record's rate does not state its supplying node "
                       "(uncertainty correlation)")
        elif stated != rmg.get("supplier"):
            out.append(f"supplier: record {stated!r} != RMG's {rmg.get('supplier')!r}")
        if supplier is not None and rmg.get("supplier") != supplier:
            out.append(f"supplier: RMG's {rmg.get('supplier')!r} != the certified walk's "
                       f"{supplier!r}")
    if rmg.get("sources") is not None:
        # absent keys compare as None (a training match has no "exact")
        for key in ("exact", "rules", "training"):
            ours, theirs = _contributors(prs.get(key)), _contributors(rmg["sources"].get(key))
            if ours != theirs:
                out.append(f"{key}: record contributors/weights {ours} != RMG {theirs}")
    return out


def _contributors(value):
    """Entry identities (index, label) and weights of a contributor list."""
    if not isinstance(value, list):
        return value

    def ident(entry):
        entry = entry or {}
        return (entry.get("index"), entry.get("label"))
    out = []
    for item in value:
        out.append((ident(item.get("rule")) if "rule" in item else None,
                    ident(item.get("entry")), item.get("reverse"),
                    float(item.get("weight")) if item.get("weight") is not None else None))
    return sorted(out, key=repr)


def _find_channel(channel_set, cert_channel):
    for ch in channel_set.channels:
        if list(ch.pathway) == cert_channel["pathway"] and \
                list(ch.products) == cert_channel["products"] and \
                ch.direction == cert_channel["direction"]:
            return ch
    return None


def _resonance_forms(molecule):
    """``molecule`` and its resonance forms (atom ids kept, nothing filtered)."""
    from rmgpy.molecule.resonance import generate_resonance_structures
    if len(molecule.atoms) < 3:
        return [molecule]
    base = molecule.copy(deep=True)
    for atom, original in zip(base.atoms, molecule.atoms):
        atom.id = original.id
    try:
        forms = generate_resonance_structures(base, keep_isomorphic=True,
                                              filter_structures=False)
    except Exception:
        forms = []
    return [molecule] + [f for f in forms
                         if sorted(a.id for a in f.atoms) == sorted(a.id for a in molecule.atoms)]


def _same_species_pinned(a, b, pinned) -> bool:
    """Whether ``a`` and ``b`` are one species (any resonance form of either,
    electrons and bond orders compared) by an isomorphism that maps every
    pinned atom id of ``a`` to its pinned id in ``b``: ``a`` against every
    form of ``b``, then every other form of ``a`` against ``b``.  Forms are
    generated only when the given structures do not already match."""
    def iso(fa, fb):
        x, y = _id_copy(fa), _id_copy(fb)
        a_by = {atom.id: atom for atom in x.atoms}
        b_by = {atom.id: atom for atom in y.atoms}
        initial = {a_by[i]: b_by[j] for i, j in pinned}
        return x.is_isomorphic(y, initial_map=initial or None, strict=True)

    if iso(a, b):
        return True
    b_forms = _resonance_forms(b)
    if any(iso(a, fb) for fb in b_forms[1:]):
        return True
    return any(iso(fa, b) for fa in _resonance_forms(a)[1:])


def _id_copy(molecule):
    out = molecule.copy(deep=True)
    for atom, original in zip(out.atoms, molecule.atoms):
        atom.id = original.id
    return out


def _move_matches(move, others, owner_of, pins) -> bool:
    """The move's products equal the record's other-side species under the
    pinned heavy-atom map (hydrogens free; any resonance form)."""
    prods = [rm.molecule_from_rows(p) for p in rm.split_rows(move.product_rows())]
    return _components_match(prods, others, owner_of, pins)


def _components_match(prods, others, owner_of, pins) -> bool:
    """Molecules ``prods`` are the species ``others`` one to one under the
    pinned atom ids ``pins`` (prods id -> others id; hydrogens free)."""
    if len(prods) != len(others):
        return False
    used = set()
    free = []
    for pm in prods:
        pinned = [(a.id, pins[a.id]) for a in pm.atoms if a.id in pins]
        if not pinned:
            free.append(pm)
            continue
        targets = {owner_of.get(b) for _a, b in pinned}
        if len(targets) != 1 or None in targets:
            return False
        (t,) = targets
        if t in used:
            return False
        used.add(t)
        if not _same_species_pinned(pm, others[t], pinned):
            return False
    for pm in free:
        hit = next((t for t in range(len(others)) if t not in used
                    and _same_species_pinned(pm, others[t], [])), None)
        if hit is None:
            return False
        used.add(hit)
    return True


def _template_fields(record):
    reactant_field = "rate_witness_reactant_graphs" if record.get(
        "rate_witness_reactant_graphs") else "reactant_graphs"
    product_field = "rate_witness_product_graphs" if record.get(
        "rate_witness_product_graphs") else "product_graphs"
    if record.get("orientation", "as_generated") == "as_generated":
        return reactant_field, product_field
    return product_field, reactant_field


def record_pins(record, side, mols, others) -> dict:
    """Template-side witness heavy atom id -> other-side witness heavy atom id
    of the record's mapped event: the record's ``atom_map`` (executable
    heavy-atom positions, reactants -> products), through the padding
    projection when the template side is a rate witness."""
    record = _as_dict(record)
    amap = {int(k): int(v) for k, v in (record.get("atom_map") or {}).items()}
    heavy_counts = [sum(1 for m in side for a in m.atoms if not a.is_hydrogen())
                    for side in (witness_molecules(record.get("reactant_graphs") or []),
                                 witness_molecules(record.get("product_graphs") or []))]
    problem = heavy_map_problem(amap, *heavy_counts)
    if problem:
        raise CertificateError(problem)
    if not amap:
        return {}   # both sides without heavy atoms (H + H -> H2)
    as_generated = record.get("orientation", "as_generated") == "as_generated"
    template_heavy = [a.id for m in mols for a in m.atoms if not a.is_hydrogen()]
    other_heavy = [a.id for m in others for a in m.atoms if not a.is_hydrogen()]
    if not side.startswith("rate_witness"):
        r_ids, p_ids = (template_heavy, other_heavy) if as_generated else \
            (other_heavy, template_heavy)
        try:
            pairs = [(r_ids[i], p_ids[j]) for i, j in amap.items()]
        except IndexError:
            raise CertificateError("the record's atom map indexes beyond its graphs") from None
    else:
        projection = (record.get("proxy_padding") or {}).get("executable_to_witness_projection")
        if not isinstance(projection, dict):
            raise CertificateError("padded record without an executable-to-witness projection")
        w_r, w_p = (mols, others) if as_generated else (others, mols)

        def side_map(name, executable_graphs, witness_mols):
            heavy_counts = [sum(1 for a in m.atoms if not a.is_hydrogen())
                            for m in witness_molecules(executable_graphs)]
            starts = [sum(heavy_counts[:k]) for k in range(len(heavy_counts))]
            w_heavy = [[a.id for a in m.atoms if not a.is_hydrogen()] for m in witness_mols]
            out = {}
            for item in projection.get(name, []):
                p = int(item["participant_index"])
                wp = int(item.get("witness_participant_index", p))
                out[starts[p] + int(item["executable_atom_index"])] = \
                    w_heavy[wp][int(item["witness_atom_index"])]
            return out
        try:
            r_map = side_map("reactants", record["reactant_graphs"], w_r)
            p_map = side_map("products", record["product_graphs"], w_p)
            pairs = [(r_map[i], p_map[j]) for i, j in amap.items()]
        except (KeyError, IndexError, ValueError, TypeError):
            raise CertificateError("the padding projection does not cover the record's atom "
                                   "map") from None
    pins = {a: b for a, b in pairs} if as_generated else {b: a for a, b in pairs}
    if len(set(pins.values())) != len(pins):
        raise CertificateError("the record's atom map is not injective")
    return pins


def heavy_map_problem(amap, reactant_heavy: int, product_heavy: int):
    """Why ``amap`` is not a complete canonical heavy-atom map (None: it is):
    a bijection of the executable reactant heavy-atom positions onto the
    product ones.  The empty map is complete only when both sides have no
    heavy atom."""
    if reactant_heavy != product_heavy:
        return (f"heavy map: {reactant_heavy} reactant and {product_heavy} product heavy atoms "
                "(not a bijection)")
    if set(amap) != set(range(reactant_heavy)) or \
            sorted(amap.values()) != list(range(product_heavy)):
        return (f"heavy map: the atom map covers {len(amap)} of {reactant_heavy} heavy atoms "
                "(a complete canonical heavy-atom bijection is required)")
    return None


def runtime_operations(record) -> list:
    """The executable operation stream the runtime applies to ``record``:
    its ``bond_ops``, then every ``feature_ops`` entry not already among them
    (``state.KMCState._record_operations``, called as is)."""
    from rmgpy.kmc.state import KMCState
    return KMCState._record_operations(None, record)


# the atom fields the runtime's set_* operations write (state.KMCState._rewrite_graph)
_RUNTIME_FIELDS = {"set_radical": "radical", "set_charge": "charge",
                   "set_lone_pairs": "lone_pairs", "set_implicit_hydrogens": "implicit_hydrogens"}


def _runtime_graph(graphs) -> dict:
    """One runtime atom graph (``state._parse_adjacency`` nodes: element,
    radical, charge, lone pairs, implicit hydrogens, edges with their
    orders) of the concatenated participants; keys are the record's atom
    positions, the indices of its operations."""
    from rmgpy.kmc.state import _parse_adjacency
    graph, offset = {}, 0
    for text in graphs:
        part = _parse_adjacency(text)
        for i, node in part.items():
            node = dict(node)
            node["edges"] = {j + offset: o for j, o in node["edges"].items()}
            graph[i + offset] = node
        offset += len(part)
    return graph


def runtime_replay(record):
    """(graph, problems): the record's reactant graphs rewritten by its
    complete executable operation stream (``runtime_operations``) with the
    runtime's semantics, an exact mirror of the graph edits of
    ``state.KMCState._rewrite_graph``: break needs the bond, form needs its
    absence, change needs the bond, form/change store ``state._order`` of the
    requested order, set_* stores the requested integer in the runtime field
    (radical, charge, lone pairs, implicit hydrogens); anything else is
    refused.  No chemistry normalisation is applied."""
    from rmgpy.kmc.state import _order
    record = _as_dict(record)
    graph = _runtime_graph(record.get("reactant_graphs") or [])
    problems = []
    for k, op in enumerate(runtime_operations(record)):
        action = op.get("action")
        try:
            if action in ("break", "form", "change"):
                first, second = (int(index) for index in op["atoms"])
                if first not in graph or second not in graph:
                    raise KeyError("an operation atom is not a reactant atom")
                existing = graph[first]["edges"].get(second)
                if action == "break":
                    if existing is None:
                        problems.append(f"operation {k}: cannot break a missing bond")
                        continue
                    graph[first]["edges"].pop(second, None)
                    graph[second]["edges"].pop(first, None)
                    continue
                if action == "form" and existing is not None:
                    problems.append(f"operation {k}: cannot form an existing bond")
                    continue
                if action == "change" and existing is None:
                    problems.append(f"operation {k}: cannot change a missing bond")
                    continue
                order = _order(op["order"])
                graph[first]["edges"][second] = order
                graph[second]["edges"][first] = order
            elif isinstance(action, str) and action.startswith("set_"):
                field = _RUNTIME_FIELDS.get(action)
                if field is None:
                    problems.append(f"operation {k}: unknown feature rewrite {action!r}")
                    continue
                atom = int(op["atom"])
                if atom not in graph:
                    raise KeyError("the operation atom is not a reactant atom")
                graph[atom][field] = int(op["value"])
            else:
                problems.append(f"operation {k}: unknown rewrite action {action!r}")
        except (KeyError, TypeError, ValueError) as exc:
            problems.append(f"operation {k} does not apply to the reactant graphs: {exc}")
    return graph, problems


def _runtime_signature(node):
    return (node.get("element"), int(node.get("radical", 0)), int(node.get("charge", 0)),
            int(node.get("lone_pairs", 0)), int(node.get("implicit_hydrogens", 0)))


def _runtime_components(graph) -> list:
    unseen, out = set(graph), []
    while unseen:
        stack, members = [min(unseen)], set()
        while stack:
            current = stack.pop()
            if current in members:
                continue
            members.add(current)
            stack.extend(set(graph[current]["edges"]) - members)
        unseen -= members
        out.append(members)
    return out


def _exact_component_match(graph, members, product, pins) -> bool:
    """``graph`` restricted to ``members`` equals runtime graph ``product``
    exactly (every node field, every edge with its order) by a bijection
    that sends each pinned atom to its pin (hydrogens free)."""
    from rmgpy.kmc.state import _order
    if len(members) != len(product):
        return False
    pinned = set(pins.values())
    candidates = {}
    for src in members:
        node = graph[src]
        options = [pins[src]] if src in pins else [t for t in product if t not in pinned]
        candidates[src] = [t for t in options if t in product
                           and _runtime_signature(product[t]) == _runtime_signature(node)
                           and len(product[t]["edges"]) == len(node["edges"])]
        if not candidates[src]:
            return False
    order = sorted(members, key=lambda a: (len(candidates[a]), a))

    def visit(position, matched, used):
        if position == len(order):
            return True
        src = order[position]
        for target in candidates[src]:
            if target in used:
                continue
            if all(_order(product[target]["edges"].get(matched[n])) == _order(o)
                   for n, o in graph[src]["edges"].items() if n in matched):
                matched[src] = target
                used.add(target)
                if visit(position + 1, matched, used):
                    return True
                del matched[src]
                used.discard(target)
        return False

    return visit(0, {}, set())


def _normalisation_problems(graph, members) -> list:
    """Requested values of a replayed component that RMG's chemistry
    normalisation (``Molecule.update``: lone pairs and charges recomputed,
    atom types, multiplicity) would change; an invalid graph is a problem."""
    from rmgpy.molecule.molecule import Atom, Bond, Molecule
    mol, atoms = Molecule(), {}
    try:
        for i in sorted(members):
            node = graph[i]
            atoms[i] = Atom(element=node["element"], radical_electrons=int(node["radical"]),
                            charge=int(node["charge"]), lone_pairs=int(node["lone_pairs"]))
            mol.add_atom(atoms[i])
        for i in sorted(members):
            for j, order in graph[i]["edges"].items():
                if i < j:
                    mol.add_bond(Bond(atoms[i], atoms[j], order=float(order)))
        normalised = mol.copy(deep=True)
        normalised.update(log_species=False, sort_atoms=False)
    except Exception as exc:  # an invalid graph is a problem, never a pass
        return [f"the replayed graph is not a valid molecule: {type(exc).__name__}: {exc}"]
    out = []
    for (i, before), after in zip(atoms.items(), normalised.atoms):
        for name in ("radical_electrons", "charge", "lone_pairs"):
            if getattr(before, name) != getattr(after, name):
                out.append(f"atom {i}: the requested {name} {getattr(before, name)} is changed "
                           f"to {getattr(after, name)} by chemistry normalisation "
                           "(Molecule.update)")
    return out


def rewrite_problems(record) -> list:
    """Replay the record's complete executable operation stream (``bond_ops``
    and ``feature_ops``, ``runtime_operations``) on its reactant graphs with
    the runtime's own semantics (``runtime_replay``; graph positions are the
    operation indices) and require, BEFORE any chemistry normalisation, the
    declared products exactly (every runtime field: element, radical
    electrons, charge, lone pairs, implicit hydrogens; every bond with its
    order) under the record's complete heavy-atom map (heavy atoms pinned,
    hydrogens free; the empty map for hydrogen-only sides), and that no
    requested value is one chemistry normalisation would change.  Empty: the
    rewrite is the record's."""
    from rmgpy.kmc.state import _parse_adjacency
    record = _as_dict(record)
    try:
        # the runtime's own reading of the graphs (no chemistry perception)
        graph, problems = runtime_replay(record)
        targets = [_parse_adjacency(text) for text in record.get("product_graphs") or []]
        amap = {int(k): int(v) for k, v in (record.get("atom_map") or {}).items()}
    except Exception as exc:  # unparsable graphs are a problem, never a pass
        return [f"the record's graphs do not parse: {type(exc).__name__}: {exc}"]
    r_heavy = [i for i in sorted(graph) if graph[i]["element"] != "H"]
    p_heavy = [(t, i) for t, part in enumerate(targets) for i in sorted(part)
               if part[i]["element"] != "H"]
    problem = heavy_map_problem(amap, len(r_heavy), len(p_heavy))
    if problem:
        return [problem]
    if not runtime_operations(record):
        return ["the record has no executable operations (bond_ops, feature_ops)"]
    if problems:
        return [f"the operations do not apply to the reactant graphs: {p}" for p in problems]
    pins = {r_heavy[i]: p_heavy[j] for i, j in amap.items()}
    used = set()
    components = _runtime_components(graph)
    if len(components) != len(targets):
        return [f"the replayed operations make {len(components)} components, the record "
                f"declares {len(targets)} products"]
    for members in sorted(components, key=lambda m: not (set(m) & set(pins))):
        owners = {pins[a][0] for a in members if a in pins}
        if len(owners) > 1:
            return ["the replayed operations join atoms the heavy-atom map sends to different "
                    "products"]
        options = list(owners) if owners else [t for t in range(len(targets)) if t not in used]
        local = {a: pins[a][1] for a in members if a in pins}
        hit = next((t for t in options if t not in used
                    and _exact_component_match(graph, members, targets[t], local)), None)
        if hit is None:
            return ["the replayed operations (bond_ops + feature_ops, runtime semantics) do not "
                    "reproduce the declared products exactly under the heavy-atom map (compared "
                    "before chemistry normalisation)"]
        used.add(hit)
        problems.extend(_normalisation_problems(graph, members))
    return problems


def binding_evidence(record) -> dict:
    """The record fields a binding verified, as stored in the binding (the
    validator requires them equal to the record's): its rewrite (bond and
    feature operations), heavy-atom map and degeneracies."""
    record = _as_dict(record)
    return {"bond_ops_sha256": _sha(record.get("bond_ops") or []),
            "feature_ops_sha256": _sha(record.get("feature_ops") or []),
            "atom_map_sha256": _sha(sorted([int(k), int(v)] for k, v in
                                           (record.get("atom_map") or {}).items())),
            "degeneracy": record.get("degeneracy"),
            "raw_path_degeneracy": record.get("raw_path_degeneracy")}


def _bound_channels(cert, binding):
    """(index, channel) of every certified channel a binding matched."""
    seen = sorted({int(m["channel"]) for m in (binding or {}).get("members", [])})
    channels = (cert or {}).get("channels", [])
    return [(i, channels[i]) for i in seen if 0 <= i < len(channels)]


def derived_certified(cert, binding) -> bool:
    """``certified`` as derived from a certificate and its binding."""
    if cert is None or (binding or {}).get("status") != "bound":
        return False
    bound = _bound_channels(cert, binding)
    return bool(cert.get("complete") and cert.get("eligible") and bound
                and not cert.get("depository_exclusion_refused")
                and all(_publishable(c) for _i, c in bound))


def _publishable(channel) -> bool:
    return channel.get("support") == chm.SUPPORTED and \
        (channel.get("eligibility") or {}).get("disposition") == "eligible"


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
                if not mp.get("mapped") or mp.get("edits") is None:
                    gaps.append(f"{name}: member {k}: a mapping without its identity "
                                "(mapped atoms, edits)")
                if mp.get("descent") is None:
                    gaps.append(f"{name}: member {k}: {mp.get('gap')}")
                    continue
                descents = cert["descents"].get(mp["descent"])
                if descents is None:
                    gaps.append(f"{name}: member {k}: descent {mp['descent']} not stored")
                    continue
                if mp["leaves"] != list(ch["pathway"]):
                    gaps.append(f"{name}: member {k}: recorded leaves {mp['leaves']}")
                for d in descents:
                    if d.get("skipped"):
                        continue
                    if d.get("root") is None:
                        gaps.append(f"{name}: descent from {d['top']} without root check")
                    table = {row[0] for row in d["atoms"]}
                    for p in iter_descent_predicates(d):
                        if p["node"] not in cert["nodes"]:
                            gaps.append(f"{name}: predicate {p['node']} without node semantics")
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
    return sorted(set(gaps), key=gaps.index)


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


def _predicate_signature(pred, path):
    out = [(path, pred["node"], pred["result"], pred.get("logic"), pred.get("invert"),
            pred.get("evaluated"), pred.get("of"), pred.get("anchoring"))]
    for k, c in enumerate(pred.get("components") or []):
        out.extend(_predicate_signature(c, path + (k,)))
    return tuple(out)


def descent_signature(descents) -> tuple:
    """The selection a mapping's descents state: every predicate tried (nested
    logic components included, in order) with its verdict, every step's
    choice and rule, the leaves.  Ids are not part of it."""
    out = []
    for d in descents:
        if d.get("skipped"):
            out.append(("skipped", d["top"], d["structure"]))
            continue
        root = _predicate_signature(d["root"], ("root",)) if d.get("root") else None
        steps = tuple((s["node"], s.get("choice"), s.get("rule"),
                       tuple(x for k, c in enumerate(s["children"])
                             for x in _predicate_signature(c, (s["node"], k))))
                      for s in d["steps"])
        out.append((d["top"], d["structure"], root, steps, d["leaf"]))
    return tuple(out)


def _descent_signature(descents, transport=None):   # backward-compatible name
    return descent_signature(descents)


def _channel_selection(ch):
    sel = ch.get("selection") or {}
    dep = ch.get("depository") or {}
    return (tuple(ch["pathway"]), ch["kind"], ch.get("source_kind"), ch.get("support"),
            canonical_json(sel.get("supplying_node")),
            canonical_json([c.get("template") for c in sel.get("contributors", [])]),
            canonical_json(_contributors((sel.get("sources") or {}).get("rules"))),
            canonical_json(_contributors((sel.get("sources") or {}).get("training"))),
            canonical_json((sel.get("sources") or {}).get("exact")),
            canonical_json(sel.get("rmg_entry")),
            canonical_json([w.get("moved") for w in sel.get("parent_walk", [])]),
            (dep.get("selected") or {}).get("index"))


def _iter_mappings(cert):
    for ch in cert["channels"]:
        for mem in ch.get("members", []):
            for mp in mem["mappings"]:
                yield ch, mem, mp


def mapping_view(cert, transport: Optional[Dict[int, int]] = None):
    """{mapping identity (``mapping_key``, ids transported): sorted tuple of
    (channel selection, descent signature)} over every certified mapping
    (``transport``: old id -> new id).  Mappings are never merged: the key
    keeps template slots, the full mapped correspondence and the edits."""
    t = _transport_fn(transport)
    out = defaultdict(list)
    for ch, _mem, mp in _iter_mappings(cert):
        sig = descent_signature(cert["descents"][mp["descent"]]) if mp.get("descent") else None
        out[mapping_key(mp, t)].append((_channel_selection(ch), sig))
    return {k: tuple(sorted(v, key=repr)) for k, v in out.items()}


def selection_view(cert, transport: Optional[Dict[int, int]] = None):
    """The whole certificate's selection, ids transported (``transport``) or
    raw: the set of (mapping identity, selections)."""
    return frozenset(mapping_view(cert, transport).items())


def _t_list(values, t):
    return tuple(sorted((t(v) for v in values), key=repr))


def _predicate_support(pred, t):
    if pred.get("logic"):
        return None
    emb = pred.get("embeddings") or {}
    guards = {k: (g.get("structure"), g.get("read")) for k, g in
              sorted((pred.get("guards") or {}).items())}
    # embeddings compared independently of the id order: their number, every
    # group atom's image set and the image union (never the raw-id ``least``)
    return (pred["node"], pred.get("anchoring"),
            tuple(sorted((lab, t(a)) for lab, a in (pred.get("anchors") or {}).items())),
            tuple(sorted((lab, _t_list(ids, t)) for lab, ids in
                         (pred.get("candidates") or {}).items())),
            tuple(sorted(((t(i), d) for i, d in pred.get("ball", [])), key=repr)),
            (emb.get("count"),
             tuple(sorted((int(g), _t_list(ids, t)) for g, ids in emb.get("by_group", []))),
             _t_list(emb.get("images", []), t)),
            canonical_json(guards), _t_list(pred.get("ignored", []), t))


def _row_support(row, t):
    return (t(row[0]),) + tuple(row[1:8]) + (
        tuple(sorted(((t(n), o) for n, o in row[8]), key=repr)),)


def support_view(cert, transport: Optional[Dict[int, int]] = None):
    """{mapping identity: transported support} per certified mapping: for
    every group predicate its anchors, candidates, ball, embeddings, guard
    values and ignored atoms, and every atom-table row (atom attributes and
    FULL bond lists) of the descents."""
    t = _transport_fn(transport)
    out = {}
    for _ch, _mem, mp in _iter_mappings(cert):
        if not mp.get("descent"):
            continue
        rows = []
        for d in cert["descents"][mp["descent"]]:
            if d.get("skipped"):
                continue
            preds = tuple(_predicate_support(p, t) for p in iter_descent_predicates(d))
            table = tuple(sorted((_row_support(r, t) for r in d.get("atoms", [])), key=repr))
            rows.append((d["top"], d["structure"], preds, table))
        out.setdefault(mapping_key(mp, t), []).append(tuple(rows))
    return {k: tuple(sorted(v, key=repr)) for k, v in out.items()}


def compare_certificates(old, new, transport) -> dict:
    """Selection and support of ``new`` against ``old`` transported by
    ``transport`` (old id -> new id): ``same_selection`` (non-empty, every
    old mapping identity present with equal selection), ``same_support``
    (equal transported support), ``fixed_point`` (both); the differences."""
    old_view = mapping_view(old, transport)
    new_view = mapping_view(new)
    missing = sorted(repr(k) for k in old_view if k not in new_view)
    changed = sorted(repr(k) for k in old_view if k in new_view and new_view[k] != old_view[k])
    # mappings of the recertified cohort (the channels the new certificate
    # selected through the old mappings' labels) that no old mapping maps to
    added = sorted(repr(k) for k in new_view if k not in old_view)
    same_selection = bool(old_view) and not missing and not changed
    old_sup, new_sup = support_view(old, transport), support_view(new)
    support_changes = sorted(repr(k)[:200] for k in old_sup
                             if k in new_view and new_sup.get(k) != old_sup[k])
    same_support = same_selection and not support_changes
    return {"same_selection": same_selection, "same_support": same_support,
            "same_cohort": not added,
            "fixed_point": same_selection and same_support and not added,
            "missing_mappings": missing, "changed_mappings": changed,
            "added_mappings": added, "support_changes": support_changes,
            "flips": flipped_predicates(old, new, transport)}


def record_certificate_refusals(record, record_certificate, rate_owner=None,
                                certifier=None) -> list:
    """Why ``record_certificate`` cannot stand for ``record`` in a
    recertification (empty: it can).  A supplied certificate undergoes the
    artifact check (``stored_certificate_problems``, one code path):
    schema, ownership against ``record`` AND the supplied rate-owning record
    (``rate_owner``: the record itself for an RMG family estimate, its
    ``reverse_of`` partner for a reference-thermo reverse; its id and rate
    source), and re-derivation against ``certifier``'s loaded database
    (none: unverified, refused); then it must be successfully bound and
    certified, with an inner certificate."""
    record = _as_dict(record)
    rc = record_certificate if isinstance(record_certificate, dict) else {}
    if rc.get("schema") != RECORD_CERTIFICATE_SCHEMA:
        return ["not a record certificate (schema): record-level recertification needs the "
                "record's owned, bound and certified record certificate"]
    by_id = {record.get("event_id"): record}
    if rate_owner is not None:
        ro = _as_dict(rate_owner)
        by_id[ro.get("event_id")] = ro
    out = [f"the record certificate: {p}"
           for p in stored_certificate_problems(record, rc, by_id, certifier)]
    binding = rc.get("binding") or {}
    if binding.get("status") != "bound":
        out.append(f"the record certificate is not bound: {binding.get('mismatches', [])[:2]}")
    if not rc.get("certified") or not derived_certified(rc.get("certificate"), binding):
        out.append(f"the record certificate is not certified: {rc.get('gaps', [])[:2]}")
    if rc.get("certificate") is None:
        out.append("the record certificate has no certificate (no RMG selection)")
    return out


def old_certificate_refusals(cert) -> list:
    """Why a certificate cannot take part in recertification (empty: it can)."""
    out = []
    if not any(True for _ in _iter_mappings(cert)):
        out.append("the certificate has no certified mapping (nothing to recertify)")
    if not cert.get("complete"):
        out.append(f"the certificate is not complete: {cert.get('gaps', [])[:2]}")
    if not cert.get("eligible"):
        why = cert.get("input_refusals") or [
            f"{';'.join(ch['pathway'])}: {ch.get('eligibility')}"
            for ch in cert.get("channels", [])
            if (ch.get("eligibility") or {}).get("disposition") != "eligible"]
        out.append(f"the certificate is not domain-eligible: {why[:2]}")
    if cert.get("depository_exclusion_refused"):
        out.append("a long-chain analogue may match a depository entry (depository exclusion "
                   f"refused for {cert['depository_exclusion_refused']})")
    return out


def transport_problems(old_cert, witness, transport) -> list:
    """Every way ``transport`` fails to carry the old witness onto
    ``witness``: not injective, an image absent, a surviving atom's element,
    isotope, radical electrons, charge or lone pairs changed, a bond among
    surviving atoms (hydrogens included) added or lost, an unmapped old atom
    that is heavy or belongs to a certified mapping.  Bond orders are not
    compared: with explicit hydrogens and the electrons compared they are a
    resonance choice, and each certificate enumerates every form.  That
    holds inside the declared v1 PS domain (``domain.check``: saturated
    backbones, phenyls, declared ends); outside it a heavy bond order can
    change with the electrons unchanged and this check does not see it."""
    problems = []
    old_atoms = {row[0]: row for row in old_cert.get("witness", {}).get("atoms", [])}
    if not old_atoms:
        return ["the old certificate does not record its witness atoms"]
    new_atoms = {row[0]: row for row in _witness_atoms([aromatic_root(m) for m in witness])}
    if len(set(transport.values())) != len(transport):
        problems.append("not injective")
    involved = set()
    for _ch, _mem, mp in _iter_mappings(old_cert):
        involved.update(a for _s, _l, a in mp["labels"])
        involved.update(a for _s, _i, _j, a in mp.get("mapped", []))
        edits = mp.get("edits") or {}
        involved.update(a for a, _b, _c in edits.get("atoms", []))
        for a, b, _o1, _o2 in edits.get("bonds", []):
            involved.update((a, b))
    for a, row in sorted(old_atoms.items()):
        if a not in transport:
            if row[1] != "H" or a in involved:
                problems.append(f"old atom {a} ({row[1]}) is not mapped")
            continue
        image = new_atoms.get(transport[a])
        if image is None:
            problems.append(f"old atom {a} maps to absent atom {transport[a]}")
            continue
        if len(row) < 8:
            problems.append(f"old atom {a}: the old witness does not record full attributes")
            continue
        changed = [name for name, k in (("element", 1), ("isotope", 4),
                                        ("radical electrons", 2), ("charge", 5),
                                        ("lone pairs", 6)) if image[k] != row[k]]
        if changed:
            problems.append(f"old atom {a} ({row[1]}, isotope {row[4]}, u{row[2]}) maps to "
                            f"{image[0]} ({image[1]}, isotope {image[4]}, u{image[2]}): "
                            f"{', '.join(changed)} changed")
            continue
        # bonds: the explicit-H constitution (bond orders are a resonance
        # choice, fixed by the constitution and the electrons compared above;
        # the aromatic root of a witness may write a ring Kekule on one side
        # and aromatic on the other)
        survivors = {transport[x] for x in old_atoms if x in transport}
        old_bonds = {transport[n] for n in row[3] if n in transport}
        new_bonds = set(image[3]) & survivors
        if old_bonds != new_bonds:
            problems.append(f"bonds of old atom {a} are not transported")
    extra = sorted(set(transport) - set(old_atoms))
    if extra:
        problems.append(f"the map names atoms {extra[:5]} outside the old witness")
    return problems


def _refused_recertification(refusals) -> dict:
    return {"certificate": None, "same_selection": False, "same_support": False,
            "same_cohort": False, "fixed_point": False, "missing_mappings": [],
            "changed_mappings": [], "added_mappings": [], "support_changes": [], "flips": [],
            "refusals": list(refusals), "reason": "; ".join(refusals)}


def flipped_predicates(old, new, transport):
    """(transported mapping labels, node path, old result, new result) for
    every certified mapping of ``old`` whose predicate verdict differs on the
    mapping of ``new`` with the same identity, or that was tried in only one."""
    def index(cert, tr):
        out = defaultdict(dict)
        for _ch, _mem, mp in _iter_mappings(cert):
            if not mp.get("descent"):
                continue
            key = mapping_key(mp, tr)
            for d in cert["descents"][mp["descent"]]:
                if d.get("skipped"):
                    continue
                for entry in descent_signature([d])[0][2] or ():
                    out[key][(d["top"], entry[0], entry[1])] = entry[2]
                for step in d["steps"]:
                    for k, c in enumerate(step["children"]):
                        for entry in _predicate_signature(c, (step["node"], k)):
                            out[key][(d["top"], entry[0], entry[1])] = entry[2]
        return out
    a = index(old, _transport_fn(transport))
    b = index(new, lambda x: x)
    flips = []
    for key in sorted(a, key=repr):
        pa, pb = a.get(key, {}), b.get(key, {})
        labels = sorted([lab, str(x)] for _s, lab, x in key[0])
        for node in sorted(set(pa) | set(pb), key=repr):
            ra, rb = pa.get(node, "not tried"), pb.get(node, "not tried")
            if ra != rb:
                flips.append({"labels": labels, "top": node[0], "path": list(map(str, node[1])),
                              "node": node[2], "old": ra, "new": rb})
    return flips


# ---------------------------------------------------------------------------
# record binding, references, certification mode
# ---------------------------------------------------------------------------
def record_owner(record) -> dict:
    """The owner block of a record certificate: the record's identity (event
    id, which hashes every semantic field), family, template, orientation,
    rate-source kind, template side, digests of its graphs, atom map and
    rewrite (``bond_ops``) and its degeneracies (rate-owner fields are added
    by ``certify_record``)."""
    record = _as_dict(record)
    side, _other = _template_fields(record)
    graphs = {f: list(record.get(f) or []) for f in (
        "reactant_graphs", "product_graphs", "rate_witness_reactant_graphs",
        "rate_witness_product_graphs")}
    amap = sorted([int(k), int(v)] for k, v in (record.get("atom_map") or {}).items())
    return {"event_id": record.get("event_id"), "family": record.get("family"),
            "template": record.get("template"), "orientation": record.get("orientation"),
            "rate_source_kind": (record.get("rate_source") or {}).get("kind"),
            "template_side": side, "graphs_sha256": _sha(graphs),
            "atom_map_sha256": _sha(amap),
            "bond_ops_sha256": _sha(record.get("bond_ops") or []),
            "degeneracy": record.get("degeneracy"),
            "raw_path_degeneracy": record.get("raw_path_degeneracy")}


def certificate_reference(record_certificate: dict) -> dict:
    """The serialized reference an EventRecord carries (``rate_certificate``):
    the sha256 of the record certificate (canonical JSON), where it lives,
    and a summary DERIVED from the stored certificate (validation recomputes
    it and requires equality)."""
    cert = record_certificate.get("certificate")
    owner = record_certificate.get("owner") or {}
    kinds = sorted({ch.get("source_kind") for ch in (cert or {}).get("channels", [])} - {None})
    return {"schema": REFERENCE_SCHEMA,
            "sha256": certificate_sha256(record_certificate),
            "store": "artifact.rate_certificates",
            "owner_event_id": owner.get("event_id"),
            "template_side": record_certificate.get("template_side"),
            "padded": bool(record_certificate.get("padded")),
            "complete": bool(record_certificate.get("complete")),
            "eligible": bool(record_certificate.get("eligible")),
            "bound": (record_certificate.get("binding") or {}).get("status") == "bound",
            "certified": bool(record_certificate.get("certified")),
            "source_kinds": kinds or [record_certificate.get("source_kind")],
            "whole_species": any(ch.get("whole_species") for ch in (cert or {}).get("channels", []))}


def certification_declaration() -> dict:
    """The artifact-level declaration of a certified compile."""
    return {"mode": CERTIFICATION_MODE, "coverage": CERTIFICATION_COVERAGE,
            "certificate_version": CERTIFICATE_VERSION,
            "record_certificate_schema": RECORD_CERTIFICATE_SCHEMA,
            "reference_schema": REFERENCE_SCHEMA,
            "required_population": "every record that fires (status enabled or irreversible) "
                                   "whose rate is an RMG family selection (an RMG family "
                                   "estimate, or a reference-thermo reverse of one) must be "
                                   "certified; a record whose certification failed is "
                                   "published refused, with its linked partner",
            "binding": "owner (event id, family, template, orientation, template side, graph, "
                       "atom-map and rewrite digests, degeneracies, rate-owning record and its "
                       "rate-source digest); binding evidence and every flag re-derived from "
                       "the stored certificate; reference summary derived from it",
            "verification": "valid only when re-derived against the loaded RMG database and "
                            "the record (family, template and tree descents resolve, the "
                            "production binding rerun equals the stored certificate); without "
                            "a database: unverified, refused"}


def requires_certification(record: dict, by_id: dict) -> bool:
    """Whether certification mode ``required`` demands that ``record`` be
    certified: it fires (status enabled or irreversible) and its rate is an
    RMG family selection (an RMG family estimate, or a reference-thermo
    reverse whose rate-owning record is one)."""
    if record.get("status", "enabled") not in FIRING_STATUSES:
        return False
    rs = record.get("rate_source") or {}
    kind = rs.get("kind")
    if kind == "RMG family estimate":
        return True
    if kind == "reference-thermo reverse":
        if rs.get("reference") == "k_library / Kc":
            return False
        owner = by_id.get(record.get("reverse_of"))
        return owner is not None and (owner.get("rate_source") or {}).get("kind") == \
            "RMG family estimate"
    return False


def owner_problems(record: dict, record_certificate: dict) -> list:
    """Why ``record_certificate``'s owner block is not ``record``'s."""
    owner = record_certificate.get("owner")
    if not isinstance(owner, dict):
        return ["the rate certificate has no owner block"]
    out = []
    for key, value in record_owner(record).items():
        if owner.get(key) != value:
            out.append(f"owner mismatch ({key}): the certificate of {owner.get('event_id')} "
                       f"is referenced by {record.get('event_id')}")
    return out


def ownership_problems(record: dict, record_certificate: dict, by_id: dict) -> list:
    """Why ``record_certificate`` does not belong to ``record`` and to its
    rate-owning record in ``by_id`` (event id -> record): the owner block,
    the rate owner's id (the record itself for an RMG family estimate, its
    ``reverse_of`` partner otherwise) and its rate source digest."""
    problems = owner_problems(record, record_certificate)
    if problems:
        return [f"rate certificate {problems[0]}"]
    owner = record_certificate["owner"]
    rate_owner = owner.get("rate_owner_event_id")
    kind = (record.get("rate_source") or {}).get("kind")
    if rate_owner is None:
        if record_certificate.get("certified"):
            return ["rate certificate owner mismatch: certified without a rate owner"]
        return []
    if kind == "RMG family estimate":
        expected_owner = record.get("event_id")
    else:
        expected_owner = record.get("reverse_of")
    if rate_owner != expected_owner or rate_owner not in by_id:
        return ["rate certificate owner mismatch (rate owner): "
                f"{rate_owner} is not the record's rate-owning record (or it was not "
                "supplied)"]
    if owner.get("rate_source_sha256") != _sha(by_id[rate_owner].get("rate_source") or {}):
        return ["rate certificate owner mismatch (rate source): the bound production "
                "selection is not the rate-owning record's"]
    return []


def check_ownership(record: dict, record_certificate: dict, by_id: dict) -> None:
    """Raise ValueError unless ``record_certificate`` belongs to ``record``."""
    problems = ownership_problems(record, record_certificate, by_id)
    if problems:
        raise ValueError(problems[0])


def _derivation_problems(record: dict, body: dict) -> list:
    """Every way the stored record certificate's claims are not derived from
    its own stored derivation (empty: consistent).  Checks the schema, the
    head against the record, the certificate's witness against the record's
    template-side graphs, the binding's coverage and evidence, and recomputes
    every flag (complete, eligible, depository refusal, publishable,
    certified)."""
    out = []
    if body.get("schema") != RECORD_CERTIFICATE_SCHEMA:
        return [f"rate certificate schema {body.get('schema')!r} is not "
                f"{RECORD_CERTIFICATE_SCHEMA!r}"]
    for key in ("event_id", "family", "template", "orientation"):
        if body.get(key) != record.get(key):
            out.append(f"rate certificate head ({key}) is not the record's")
    if body.get("rate_source_kind", (record.get("rate_source") or {}).get("kind")) != \
            (record.get("rate_source") or {}).get("kind"):
        out.append("rate certificate head (rate_source_kind) is not the record's")
    cert = body.get("certificate")
    binding = body.get("binding")
    if not isinstance(binding, dict) or not isinstance(binding.get("mismatches"), list) or \
            binding.get("status") not in ("bound", "refused", "not_applicable"):
        return out + ["rate certificate binding is malformed"]
    status, mismatches = binding["status"], binding["mismatches"]
    if status == "bound" and (mismatches or not binding.get("members")):
        out.append("rate certificate binding claims bound without members or with mismatches")
    if status == "refused" and not mismatches:
        out.append("rate certificate binding is refused without a reason")
    if status == "not_applicable" and (cert is not None or binding.get("members")):
        out.append("rate certificate binding is not applicable but binds a certificate")
    if cert is not None:
        if not isinstance(cert, dict) or cert.get("version") != CERTIFICATE_VERSION:
            return out + [f"rate certificate schema: inner certificate version "
                          f"{(cert or {}).get('version')!r} is not {CERTIFICATE_VERSION!r}"]
        if cert.get("family") != record.get("family"):
            out.append("rate certificate witness/family is not the record's")
        side = body.get("template_side")
        if side not in _template_fields(record)[:1]:
            out.append("rate certificate template side is not the record's")
        else:
            try:
                expected = _constitution(_witness_atoms(witness_molecules(record[side])))
            except Exception as exc:  # unparsable graphs are a mismatch, never a pass
                expected = f"unparsable: {exc}"
            stored = _constitution((cert.get("witness") or {}).get("atoms", []))
            if stored != expected:
                out.append("rate certificate witness is not the record's template-side graph")
        channels = cert.get("channels", [])
        for m in binding.get("members", []):
            ci, k = int(m.get("channel", -1)), int(m.get("member", -1))
            if not (0 <= ci < len(channels)) or not (
                    0 <= k < len(channels[ci].get("members", []))) or any(
                    not (0 <= int(i) < len(channels[ci]["members"][k]["mappings"]))
                    for i in m.get("mappings", [])) or not m.get("mappings"):
                out.append("rate certificate binding names a member the certificate lacks")
                break
        gaps = list(cert.get("gaps", []))
        derived_gaps = completeness_gaps(cert)
        if not channels:
            derived_gaps.append("no channel on the witness matches the selector")
        if not set(derived_gaps) <= set(gaps) or any(g not in gaps for g in derived_gaps):
            out.append("rate certificate gaps are not derived from its certificate")
        complete = not gaps
        eligible = bool(channels) and not cert.get("input_refusals") and all(
            (ch.get("eligibility") or {}).get("disposition") == "eligible" for ch in channels)
        refused = [";".join(ch["pathway"]) for ch in channels
                   if not ch.get("whole_species")
                   and exclusion_refused(ch.get("depository_exclusion"))]
        for key, value in (("complete", complete), ("eligible", eligible),
                           ("depository_exclusion_refused", refused)):
            if cert.get(key) != value:
                out.append(f"rate certificate flag {key} is not derived from its certificate")
        for ch in channels:
            if "publishable" in ch and ch["publishable"] != _publishable(ch):
                out.append("rate certificate flag publishable is not derived from its channel")
                break
    else:
        complete = eligible = False
        if status == "bound":
            out.append("rate certificate binding claims bound without a certificate")
    for key, value in (("complete", complete), ("eligible", eligible),
                       ("certified", derived_certified(cert, binding))):
        if body.get(key) != value:
            out.append(f"rate certificate flag {key} is not derived from its stored derivation")
    if status == "bound" and binding.get("evidence") != binding_evidence(record):
        diff = sorted(k for k in set(binding_evidence(record)) | set(binding.get("evidence") or {})
                      if (binding.get("evidence") or {}).get(k) != binding_evidence(record).get(k))
        out.append(f"rate certificate binding evidence ({', '.join(diff)}) is not the record's "
                   "(rewrite or degeneracy)")
    return out


def _entry_semantics(entry) -> dict:
    """The semantics the certifier records for a tree node (``_Recorder``)."""
    nodes = {}
    _Recorder(None, None, None, nodes).node_semantics(entry)
    return nodes[entry.label]


def tree_problems(family, cert, template=()) -> list:
    """Every way ``cert``'s stored derivation is not a derivation in the
    loaded family (empty: each stored claim resolves there): every label of
    ``template`` (the record's) and of every channel pathway is a node of the
    family's tree; every certified mapping's
    descent is stored and non-empty; every descent starts at a top node of
    the tree and is a real root-to-selected-node chain of it (each step's
    children are the node's children in the tree's order, a step continues
    at the previous step's choice, a matched root has steps and the leaf is
    where ``descend_tree`` stops); its leaves are the mapping's; and every
    node it evaluates exists, with the stored semantics equal to the
    semantics of the loaded node."""
    groups = getattr(family, "groups", None)
    entries = getattr(groups, "entries", None)
    if not isinstance(entries, dict):
        return [f"family {getattr(family, 'label', family)!r} has no tree in the loaded "
                "database"]
    tops = {e.label for e in groups.top}
    nodes = cert.get("nodes") or {}
    descents = cert.get("descents") or {}
    out = []

    def node_problem(label):
        entry = entries.get(label)
        if entry is None:
            return f"node {label!r} does not exist in the family tree of {family.label}"
        if nodes.get(label) != _entry_semantics(entry):
            return f"node {label!r}: the stored semantics do not resolve in the family tree"
        return None

    for ch in [{"pathway": list(template)}] + list(cert.get("channels", [])):
        for label in ch.get("pathway", []):
            if label not in entries:
                out.append(f"template label {label!r} does not exist in family {family.label} "
                           "(loaded RMG database)")
        for mem in ch.get("members", []):
            for mp in mem.get("mappings", []):
                ref = mp.get("descent")
                if ref is None:
                    continue
                chain = descents.get(ref)
                if not isinstance(chain, list) or not chain:
                    out.append(f"descent {ref!r} is empty or missing: not a root-to-node chain "
                               "of the family tree")
                    continue
                if mp.get("leaves") != descent_leaves(chain):
                    out.append(f"descent {ref!r}: the mapping's leaves are not its descent's "
                               "(family tree)")
    for label, semantics in sorted(nodes.items()):
        # tree nodes by label; anonymous logic components are keyed by their
        # repr and re-derived with their parents
        if label in entries and semantics != _entry_semantics(entries[label]):
            out.append(f"node {label!r}: the stored semantics do not resolve in the family "
                       "tree")
    for ref, chain in sorted(descents.items()):
        for d in chain if isinstance(chain, list) else []:
            top = d.get("top")
            if top not in tops:
                out.append(f"descent {ref!r} starts at {top!r}, not a top node of the family "
                           "tree")
                continue
            if d.get("skipped"):
                continue
            root = d.get("root") or {}
            steps = d.get("steps") or []
            problem = node_problem(top)
            if problem:
                out.append(problem)
            if root.get("node") != top:
                out.append(f"descent {ref!r}: the root check is not its top node's (family "
                           "tree)")
                continue
            if root.get("result") != "matched":
                if steps or d.get("leaf") is not None:
                    out.append(f"descent {ref!r}: a failed root check descends (family tree)")
                continue
            if not steps:
                out.append(f"descent {ref!r}: a matched root without steps (family tree)")
                continue
            current = top
            for k, step in enumerate(steps):
                entry = entries.get(step.get("node"))
                if step.get("node") != current or entry is None:
                    out.append(f"descent {ref!r} step {k}: {step.get('node')!r} is not the "
                               f"node the chain reached ({current!r}) in the family tree")
                    break
                children = [c.get("node") for c in step.get("children", [])]
                if children != [c.label for c in entry.children]:
                    out.append(f"descent {ref!r} step {k}: the children of {current!r} are "
                               "not the family tree's")
                    break
                for label in children:
                    problem = node_problem(label)
                    if problem:
                        out.append(problem)
                choice = step.get("choice")
                last = k == len(steps) - 1
                if choice is not None and choice not in children:
                    out.append(f"descent {ref!r} step {k}: choice {choice!r} is not a child of "
                               f"{current!r} in the family tree")
                    break
                if not last:
                    if choice is None:
                        out.append(f"descent {ref!r} step {k}: the chain continues after "
                                   "stopping (family tree)")
                        break
                    current = choice
                    continue
                expected_leaf = choice if choice is not None else current
                if choice is not None and not choice.startswith("Others-") or \
                        d.get("leaf") != expected_leaf:
                    out.append(f"descent {ref!r}: the leaf {d.get('leaf')!r} is not where the "
                               f"family tree descent stops ({expected_leaf!r})")
    return sorted(set(out), key=out.index)


def rederived_certificate(record: dict, by_id: dict, certifier: "Certifier") -> dict:
    """The record certificate re-derived from the record and the loaded
    database: ``certify_record`` (the production binding path, RMG's own tree
    evaluation and get_kinetics; nothing re-estimated) through the record's
    rate-owning ``reverse_of`` partner in ``by_id``; an exception gives the
    owned error certificate, as at compile time."""
    try:
        return certifier.certify_record(record, rate_owner=_rate_owner(record, by_id))
    except Exception as exc:  # the same verdict certify_records stores
        return _error_certificate(record, exc)


def stored_certificate_problems(record: dict, body: dict, by_id: dict,
                                certifier: Optional["Certifier"]) -> list:
    """Every reason the stored record certificate ``body`` is not valid for
    ``record`` (empty: valid).  Shared by the artifact check
    (``validate_certification``) and recertification
    (``record_certificate_refusals``).  Stored assertions are never turned
    into a verdict: after the schema and the ownership (record and rate
    owner in ``by_id``), the derivation is re-established against the loaded
    database (``certifier``): the family exists, every template label and
    every descent resolve in its tree (``tree_problems``), the stored claims
    are consistent with each other (``_derivation_problems``) and the whole
    certificate (selection, descents, node semantics, binding members, heavy
    map and rewrite evidence, flags) equals the one re-derived from the
    record (``rederived_certificate``; the owner block is the ownership
    check's).  Without a database the verdict is ``unverified``, a
    problem."""
    if not isinstance(body, dict) or body.get("schema") != RECORD_CERTIFICATE_SCHEMA:
        return [f"rate certificate schema {(body or {}).get('schema')!r} is not "
                f"{RECORD_CERTIFICATE_SCHEMA!r}"]
    problems = ownership_problems(record, body, by_id)
    if problems:
        return problems
    if certifier is None:
        return ["rate certificate unverified: no RMG database is loaded to re-derive it (a "
                "stored certificate is never valid on its own assertions)"]
    cert = body.get("certificate")
    if cert is not None:
        family = certifier.families.get(record.get("family")) \
            if hasattr(certifier.families, "get") else None
        if family is None:
            return [f"rate certificate family {record.get('family')!r} does not exist in the "
                    "loaded RMG database"]
        if not isinstance(cert, dict):
            return ["rate certificate: the inner certificate is malformed"]
        labels = [lab for lab in str(record.get("template") or "").split(";") if lab]
        try:
            problems = tree_problems(family, cert, labels)
        except Exception as exc:  # a malformed derivation is a problem, never a pass
            problems = [f"the stored descents are malformed ({type(exc).__name__}: {exc}): "
                        "not a chain of the family tree"]
        if problems:
            return [f"rate certificate {p}" for p in problems]
    problems = _derivation_problems(record, body)
    if problems:
        return problems
    fresh = rederived_certificate(record, by_id, certifier)
    stored = {k: v for k, v in body.items() if k != "owner"}
    derived = {k: v for k, v in fresh.items() if k != "owner"}
    derived["event_id"] = record.get("event_id")
    if canonical_json(stored) != canonical_json(derived):
        keys = sorted(k for k in set(stored) | set(derived)
                      if canonical_json(stored.get(k)) != canonical_json(derived.get(k)))
        return [f"rate certificate is not re-derived: the production binding rerun on the "
                f"record against the loaded database differs in {keys}"]
    return []


def validate_certification(artifact: dict, families=None) -> None:
    """Artifact-level certification check (``compiler.validate_artifact``).

    Without a ``rate_certification`` declaration the artifact must carry no
    certificate store and no reference.  With it (mode ``required``,
    coverage ``every record``): every record carries a reference that
    resolves, whose stored certificate hashes to it (digest) and is VALID
    for the record (``stored_certificate_problems``: ownership, then the
    derivation re-established against the loaded database ``families`` and
    the record itself, never read off the stored assertions), whose summary
    equals the one derived from the stored certificate; every stored
    certificate is referenced; and every record of the required population
    (``requires_certification``) is certified.  Without a loaded database
    (``families`` None) the verdict is ``unverified``, which refuses
    publication in required mode."""
    records = artifact.get("records", [])
    store = artifact.get("rate_certificates")
    declaration = artifact.get("rate_certification")
    references = [r.get("rate_certificate") for r in records]
    if declaration is None:
        if store is not None or any(ref is not None for ref in references):
            raise ValueError("rate certificates without a declared certification mode")
        return
    if declaration != certification_declaration():
        raise ValueError("unknown rate certification mode declaration")
    if not isinstance(store, dict):
        raise ValueError("certification mode required but the artifact has no certificate store")
    if families is None:
        raise ValueError("rate certification unverified: no RMG database is loaded to "
                         "re-derive the certificates; mode required refuses publication")
    certifier = Certifier(families)
    by_id = {r.get("event_id"): r for r in records}
    referenced = set()
    for record in records:
        reference = record.get("rate_certificate")
        if reference is None:
            raise ValueError(f"certification coverage: record {record.get('event_id')} has no "
                             "rate certificate (mode required: every record)")
        sha = reference.get("sha256")
        if sha not in store:
            raise ValueError("record rate certificate reference does not resolve")
        body = store[sha]
        if certificate_sha256(body) != sha:
            raise ValueError("record rate certificate digest mismatch")
        problems = stored_certificate_problems(record, body, by_id, certifier)
        if problems:
            raise ValueError(problems[0])
        if certificate_reference(body) != reference:
            raise ValueError("rate certificate reference summary contradicts the stored "
                             "certificate")
        if requires_certification(record, by_id) and not body.get("certified"):
            raise ValueError(f"certification required: record {record.get('event_id')} fires "
                             f"({record.get('status')}) with an uncertified RMG selection "
                             f"({(body.get('gaps') or ['no reason'])[0]})")
        referenced.add(sha)
    unreferenced = set(store) - referenced
    if unreferenced:
        raise ValueError(f"{len(unreferenced)} unreferenced stored rate certificate(s)")


def _rate_owner(data, by_id):
    if (data.get("rate_source") or {}).get("kind") == "reference-thermo reverse":
        return by_id.get(data.get("reverse_of"))
    return None


def _error_certificate(data, exc) -> dict:
    owner = record_owner(data)
    owner.update({"rate_owner_event_id": None, "rate_source_sha256": None})
    return {"schema": RECORD_CERTIFICATE_SCHEMA, "owner": owner,
            "event_id": data.get("event_id"), "family": data.get("family"),
            "template": data.get("template"),
            "rate_source_kind": (data.get("rate_source") or {}).get("kind"),
            "orientation": data.get("orientation"), "certificate": None,
            "complete": False, "eligible": False, "certified": False,
            "binding": {"status": "refused", "members": [], "mismatches": [
                f"{type(exc).__name__}: {exc}"]},
            "gaps": [f"{type(exc).__name__}: {exc}"], "source_kind": "error"}


def certify_records(records, families, certifier: Optional[Certifier] = None, *,
                    required: bool = True):
    """Certify every record from its actual rate source (the rate witness
    graphs when padded; a reverse record through its rate-owning
    ``reverse_of`` partner) and attach the reference.  Every record gets a
    certificate (coverage: every record); a failure is stored as an owned
    ``error`` certificate, never dropped.

    With ``required`` (the certified compile) a record of the required
    population (``requires_certification``) whose certification failed is
    published REFUSED by status ("rate certification refused: <reason>"),
    together with its linked ``reverse_of`` partner (one rate); event ids
    are recomputed and every link remapped, and each owner block is bound to
    the final record.  Returns (records, store, replacements old id -> new
    id)."""
    from dataclasses import replace
    certifier = certifier or Certifier(families)
    records = list(records)
    data = [_as_dict(r) for r in records]
    by_id = {d.get("event_id"): r for d, r in zip(data, records)}
    rcs = []
    for record, d in zip(records, data):
        try:
            rcs.append(certifier.certify_record(record, rate_owner=_rate_owner(d, by_id)))
        except Exception as exc:  # reported, never silently dropped
            rcs.append(_error_certificate(d, exc))
    replacements = {}
    if required:
        data_by_id = {d.get("event_id"): d for d in data}
        reasons = {}
        for d, rc in zip(data, rcs):
            if requires_certification(d, data_by_id) and not rc.get("certified"):
                reasons[d["event_id"]] = (rc.get("gaps") or ["not certified"])[0]
        for event_id in list(reasons):
            partner = data_by_id[event_id].get("reverse_of")
            if partner in data_by_id and partner not in reasons and \
                    data_by_id[partner].get("status", "enabled") in FIRING_STATUSES:
                reasons[partner] = f"its linked partner {event_id} failed: {reasons[event_id]}"
        if reasons:
            from rmgpy.kmc.compiler import _replace_event_ids
            from rmgpy.kmc.event_record import EventRecord
            changed = []
            for record, d in zip(records, data):
                if d["event_id"] in reasons:
                    record = replace(record, status="refused", status_reason=(
                        f"rate certification refused: {reasons[d['event_id']]}")[:2000],
                        event_id="")
                    replacements[d["event_id"]] = record.event_id
                changed.append(record)
            records = []
            for record in changed:
                relinked = _replace_event_ids(
                    {k: v for k, v in record.to_dict().items() if k != "rate_certificate"},
                    replacements)
                new_id = relinked.pop("event_id")
                relinked["event_id"] = ""
                fresh = EventRecord.from_dict(relinked)
                if fresh.event_id != new_id:
                    raise ValueError("event id changed while remapping certification links")
                records.append(fresh)
    store, out = {}, []
    final_by_id = {_as_dict(r).get("event_id"): r for r in records}
    for record, rc in zip(records, rcs):
        d = _as_dict(record)
        rc = dict(rc)
        owner = record_owner(d)
        rate_owner = (rc.get("owner") or {}).get("rate_owner_event_id")
        owner.update({"rate_owner_event_id": replacements.get(rate_owner, rate_owner),
                      "rate_source_sha256": (rc.get("owner") or {}).get("rate_source_sha256")})
        rc["owner"] = owner
        rc["event_id"] = d.get("event_id")
        if owner["rate_owner_event_id"] is not None and \
                owner["rate_owner_event_id"] not in final_by_id:
            raise ValueError("certification lost a rate-owning record")
        ref = certificate_reference(rc)
        store[ref["sha256"]] = rc
        out.append(replace(record, rate_certificate=ref) if not isinstance(record, dict)
                   else dict(record, rate_certificate=ref))
    return out, store, replacements


def attach_rate_certificates(records, families, certifier: Optional[Certifier] = None):
    """``certify_records`` without the certification gate: every record gets
    its certificate and reference, statuses and event ids unchanged.
    Returns (records, store)."""
    records, store, _replacements = certify_records(records, families, certifier,
                                                    required=False)
    return records, store


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


def certify_record(record, families, certifier: Optional[Certifier] = None,
                   rate_owner=None) -> dict:
    return (certifier or Certifier(families)).certify_record(record, rate_owner=rate_owner)


def recertify(record, witness, atom_map, *, families, certificate=None,
              certifier: Optional[Certifier] = None, rate_owner=None, domain=None) -> dict:
    return (certifier or Certifier(families)).recertify(record, witness, atom_map, certificate,
                                                        rate_owner=rate_owner, domain=domain)
