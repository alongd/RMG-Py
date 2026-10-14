"""S2 comparison harness for rmgMovesTest.py (not collected).

RMG side: ``rmg_moves.generate_moves`` over every component (unimolecular
families) or every unordered pair of distinct components (bimolecular), with
the case's reacting-atom pattern applied as ``select`` (refused mappings are
returned as ``excluded``; nothing is hidden).  The pattern is stated on the
RMG labels and evaluated on the reactant graph only (element, degree, ring
membership by RMG ring perception, heavy-atom distance), so it is the same
for every resonance form.

Matcher side: the S2 hand-written template set (``windowMatcherTest``) matched
on the same state; each transition's product rows go through the SAME RMG
converter (``rmg_moves.product_key``).

Both sets are compared by mapped product modulo resonance; every one-sided
move is listed with its atom ids and classified.
"""

import time
from collections import Counter, defaultdict
from pathlib import Path

from rmgpy.kmc import rmg_moves as rm

DATABASE = Path("/home/alon/runs/phase2b/database")
FAMILY_NAMES = ["H_Abstraction", "R_Addition_MultipleBond", "R_Recombination",
                "Disproportionation", "intra_H_migration"]
_DB = {}


def database(path=None):
    path = str(path or DATABASE)
    if path not in _DB:
        from rmgpy.data.rmg import RMGDatabase
        db = RMGDatabase()
        db.load_kinetics(path + "/input/kinetics", reaction_libraries=[], seed_mechanisms=None,
                         kinetics_families=FAMILY_NAMES, kinetics_depositories=[])
        _DB[path] = db
    return _DB[path]


def families():
    return database().kinetics.families


# --------------------------------------------------------------------------
# graph facts (RMG ring perception), form independent
# --------------------------------------------------------------------------
class Graph:
    def __init__(self, molecules):
        self.element, self.adj, self.ring, self.six_ring_c, self.mol_of = {}, {}, set(), set(), {}
        for k, mol in enumerate(molecules):
            for atom in mol.atoms:
                self.element[atom.id] = atom.element.symbol
                self.adj[atom.id] = tuple(sorted(o.id for o in atom.bonds))
                self.mol_of[atom.id] = k
                if mol.is_atom_in_cycle(atom):
                    self.ring.add(atom.id)
            for ring in mol.get_smallest_set_of_smallest_rings():
                if len(ring) == 6 and all(a.element.symbol == "C" for a in ring):
                    self.six_ring_c.update(a.id for a in ring)
        self._dist = {}

    def carbon(self, a):
        return self.element[a] == "C"

    def sp3_carbon(self, a):
        return self.carbon(a) and len(self.adj[a]) == 4

    def non_ring_heavy(self, a):
        return [b for b in self.adj[a] if self.element[b] != "H" and b not in self.ring]

    def benzylic(self, a):
        return self.carbon(a) and a not in self.ring and any(b in self.six_ring_c for b in self.adj[a])

    def chain_end(self, a):
        return self.carbon(a) and a not in self.ring and len(self.non_ring_heavy(a)) == 1

    def heavy_distance(self, a, b):
        if a not in self._dist:
            dist, frontier = {a: 0}, [a]
            while frontier:
                nxt = []
                for v in frontier:
                    for w in self.adj[v]:
                        if w not in dist and self.element[w] != "H":
                            dist[w] = dist[v] + 1
                            nxt.append(w)
                frontier = nxt
            self._dist[a] = dist
        return self._dist[a].get(b)


# --------------------------------------------------------------------------
# cases: (RMG family, direction, molecularity, core labels, scope on labels)
# --------------------------------------------------------------------------
def _lab(labels, name, slot=None):
    hits = [a for s, lab, a in labels if lab == name and (slot is None or s == slot)]
    assert len(hits) == 1, (name, labels)
    return hits[0]


def scope_ha(g, labels):
    x, y = _lab(labels, "*3"), _lab(labels, "*1")
    if not g.benzylic(x):
        return "acceptor *3 not a benzylic carbon"
    if not g.sp3_carbon(y):
        return "donor *1 not an sp3 carbon"
    return None


def scope_bs(g, labels):
    x, alpha, beta = _lab(labels, "*2"), _lab(labels, "*1"), _lab(labels, "*3")
    if not g.carbon(x) or x in g.ring:
        return "radical *2 not a non-ring carbon"
    if len([b for b in g.non_ring_heavy(x) if g.carbon(b)]) < 2:
        return "radical *2 has < 2 non-ring carbon neighbours"
    if not g.sp3_carbon(alpha) or alpha in g.ring:
        return "alpha *1 not a non-ring sp3 carbon"
    if not g.sp3_carbon(beta) or beta in g.ring:
        return "leaving *3 not a non-ring sp3 carbon"
    return None


def scope_rc(g, labels):
    for s in (0, 1):
        if not g.chain_end(_lab(labels, "*", s)):
            return f"slot {s} '*' not a chain-end carbon"
    return None


def scope_ds(g, labels):
    x, y, z, h = (_lab(labels, k) for k in ("*1", "*2", "*3", "*4"))
    if not g.benzylic(x):
        return "acceptor *1 not a benzylic carbon"
    if not g.carbon(z) or z in g.ring:
        return "donor radical *3 not a non-ring carbon"
    if not g.sp3_carbon(y) or y in g.ring:
        return "H-bearing *2 not a non-ring sp3 carbon"
    if g.element[h] != "H":
        return "*4 not H"
    return None


def scope_hs(n):
    def scope(g, labels):
        x, y = _lab(labels, "*1"), _lab(labels, "*2")
        if not g.benzylic(x):
            return "radical *1 not a benzylic carbon"
        if not g.sp3_carbon(y):
            return "donor *2 not an sp3 carbon"
        if g.heavy_distance(x, y) != n - 1:
            return f"shortest heavy path *1-*2 is {g.heavy_distance(x, y)}, not {n - 1}"
        return None
    return scope


def scope_hh(g, labels):
    for s in (0, 1):
        a = _lab(labels, "*", s)
        if g.element[a] != "H" or g.adj[a]:
            return f"slot {s} '*' not a free H atom"
    return None


# matcher family -> (RMG family, direction, molecularity, core labels, scope)
CASES = {
    "h_abstraction_benzylic": ("H_Abstraction", "forward", 2, ("*1", "*2", "*3"), scope_ha),
    "beta_scission_midchain": ("R_Addition_MultipleBond", "reverse", 1, ("*1", "*2", "*3"), scope_bs),
    "recombination_chain_end": ("R_Recombination", "forward", 2, ("*",), scope_rc),
    "disproportionation": ("Disproportionation", "forward", 2, ("*1", "*2", "*3", "*4"), scope_ds),
    "h_shift": ("intra_H_migration", "forward", 1, ("*1", "*2", "*3"), scope_hs(5)),
    "h_atom_recombination": ("R_Recombination", "forward", 2, ("*",), scope_hh),
    # the contract's set has the 1,5 shift; 1,4 and 1,6 are the neighbouring S2 classes
    "h_shift4": ("intra_H_migration", "forward", 1, ("*1", "*2", "*3"), scope_hs(4)),
    "h_shift6": ("intra_H_migration", "forward", 1, ("*1", "*2", "*3"), scope_hs(6)),
}
CONTRACT_CASES = ("h_abstraction_benzylic", "beta_scission_midchain", "recombination_chain_end",
                  "disproportionation", "h_shift", "h_atom_recombination")


def core_atoms(labels, core):
    return tuple(sorted(a for _s, lab, a in labels if lab in core))


# --------------------------------------------------------------------------
# the two sides
# --------------------------------------------------------------------------
def rmg_side(state, case, restrict=True, keys=None):
    """All RMG moves of the case on the state (dict with moves, scope, timing)."""
    fam_name, direction, arity, core, scope = CASES[case]
    family = families()[fam_name]
    mols = rm.molecules_from_state(state)
    g = Graph(mols)
    select = (lambda c: scope(g, c.labels) is None) if restrict else None
    t0 = time.perf_counter()
    moves = []
    if arity == 1:
        for mol in mols:
            if mol.multiplicity > 1:
                moves.extend(rm.generate_moves(family, [mol], direction, select=select))
    else:
        for i in range(len(mols)):
            for j in range(i + 1, len(mols)):
                if mols[i].multiplicity == 1 and mols[j].multiplicity == 1:
                    continue  # no family of the S2 set reacts two closed shells
                moves.extend(rm.generate_moves(family, [mols[i], mols[j]], direction,
                                               select=select))
    t1 = time.perf_counter()
    keys = keys or rm.product_key
    kept, filtered = defaultdict(list), defaultdict(list)
    out_of_scope = []
    for m in moves:
        reason = scope(g, m.labels)
        if reason is not None:
            out_of_scope.append((m, reason))
            continue
        if m.status == "kept":
            kept[keys(m)].append(m)
        elif m.products is not None:
            filtered[keys(m)].append(m)
    t2 = time.perf_counter()
    return {"moves": moves, "kept": dict(kept), "filtered": dict(filtered),
            "out_of_scope": out_of_scope, "graph": g,
            "generate_s": t1 - t0, "key_s": t2 - t1}


def matcher_side(state, case, keys=None):
    import windowMatcherTest as wmt
    from rmgpy.kmc import window_matcher as wm
    family, n = case, None
    if case.startswith("h_shift"):
        family, n = "h_shift", int(case[len("h_shift"):] or 5)
    templates = wmt.templates_for(family, n)
    t0 = time.perf_counter()
    result = wm.match(templates, state)
    t1 = time.perf_counter()
    keys = keys or rm.product_key
    out = defaultdict(list)
    for tr in result.transitions:
        out[keys(wm.product_rows(state, tr))].append(tr)
    t2 = time.perf_counter()
    return {"transitions": result.transitions, "by_key": dict(out), "templates": len(templates),
            "match_s": t1 - t0, "key_s": t2 - t1}


def edited_atoms(tr):
    atoms = set()
    for e in tr.edits:
        if hasattr(e, "a"):
            atoms.update((e.a, e.b))
        else:
            atoms.add(e.atom)
    return tuple(sorted(atoms))


def compare(state, case, label, restrict=True):
    """Per-case table row plus the classified one-sided moves.

    Classes: ``RMG filter`` (the matcher's product is the product of an
    in-scope RMG move with a non-kept status, and of no kept move);
    ``template scope`` (RMG keeps the product only through a mapping the
    case's scope excludes, or vice versa); ``genuine miss`` otherwise."""
    fam_name, direction, arity, core, _scope = CASES[case]
    rmg = rmg_side(state, case, restrict=restrict)
    mat = matcher_side(state, case)
    rk, mk, fk = set(rmg["kept"]), set(mat["by_key"]), set(rmg["filtered"])
    mismatches = []
    for key in sorted(mk - rk):
        tr = mat["by_key"][key][0]
        if key in fk:
            hits = rmg["filtered"][key]
            cls = "RMG filter: " + "; ".join(sorted({f"{m.status} ({m.reason})" for m in hits}))
        else:
            cls = "genuine miss (RMG): no in-scope RMG move makes this product"
        mismatches.append({"side": "matcher-only", "atoms": edited_atoms(tr),
                           "template": tr.template, "class": cls})
    for key in sorted(rk - mk):
        m = rmg["kept"][key][0]
        mismatches.append({"side": "rmg-only", "atoms": core_atoms(m.labels, core),
                           "labels": m.labels, "template": m.template,
                           "class": "genuine miss (matcher): no transition makes this product"})
    built = Counter(m.status for m in rmg["moves"])
    in_scope = Counter(m.status for m in rmg["moves"]) - Counter(
        m.status for m, _r in rmg["out_of_scope"])
    row = {
        "state": label, "case": case, "rmg_family": f"{fam_name} ({direction})",
        "rmg_all": len(rmg["moves"]), "rmg_status": dict(sorted(built.items())),
        "rmg_in_scope_status": dict(sorted(in_scope.items())),
        "rmg_kept_in_scope_moves": sum(len(v) for v in rmg["kept"].values()),
        "rmg_kept_in_scope_distinct": len(rk),
        "rmg_filtered_in_scope_distinct": len(fk - rk),
        "matcher_transitions": len(mat["transitions"]), "matcher_distinct": len(mk),
        "matched": len(rk & mk), "matcher_only": len(mk - rk), "rmg_only": len(rk - mk),
        "matcher_only_rmg_filter": sum(1 for x in mismatches if x["class"].startswith("RMG filter")),
        "genuine_misses": sum(1 for x in mismatches if x["class"].startswith("genuine")),
        "rmg_generate_s": round(rmg["generate_s"], 1), "rmg_key_s": round(rmg["key_s"], 1),
        "matcher_s": round(mat["match_s"] + mat["key_s"], 1),
    }
    scope_reasons = Counter(r for _m, r in rmg["out_of_scope"])
    row["out_of_scope_by_reason"] = dict(sorted(scope_reasons.items()))
    return row, mismatches, rmg, mat
