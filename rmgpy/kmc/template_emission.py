"""Executable local window templates from recertified closed windows (I-138).

For every compiled record whose members ``windows.py`` closed at a
recertified fixed point, this module cuts one executable template per
distinct context (``windows.run_member`` context representatives: the
all-artificial base and every context that merged with no earlier one),
weights it with the record's certified rate over the whole channel, proves
that weight local over extension vectors, and serializes it on the record
(``executable_template``, ``rewrite_record/0.2``).  It consumes windows only
through their public results (closed window and witness of each context,
the merge pairs, refusal statuses, ladder items): nothing here depends on how
the window roster is generated.  Every failure is a refusal by name or a
compile error (``EmissionError``); nothing falls back.

Template from window (item 1)
-----------------------------
``cut_template`` cuts the closed window of one context as a marked graph:

* sides: the template direction's reactant and product rows of the
  context's fixed-point witness (``reverse``: the record runs opposite to the
  member's direction, so the witness products are the template reactants),
  each written in the matcher's first resonance form
  (``window_matcher.component_forms``; the closure holds every pi system
  touching the window whole, so the edits between the two first forms are
  window edits);
* atoms: every window atom; ``closed`` when every neighbour (both sides) is
  in the window (full degree), else ``port`` with its ``exterior`` multiset
  (bond order, element, ring membership of the outside neighbour).  The
  outside neighbours and their bond orders must be equal on both sides (an
  edit never crosses the cut) and no edit may touch a port
  (``EmissionError``);
* port labels come from the witness's construction metadata
  (``boundaries.BoundaryMetadata``), never from symmetry: a port is
  ``artificial`` when its exterior reaches an artificial witness port (the
  chain continues and the fixed point proved the window does not read it),
  ``physical`` when it reaches only physical witness ports (a declared end
  just outside the window), ``interior`` when it reaches no witness port (a
  pendant outside the window).  A window atom that is itself a witness port
  carries that port (``witness_port``: participant, port id, tag,
  physicality): a physical end inside the window is a closed atom;
* participants: the template-direction reactant components the window
  touches (one or two); participant numbering, atom order and the hash are
  canonical (individualization-refinement over the marked graph, minimum
  serialization; interchangeable hydrogen leaves are not branched);
* predicates: none.  ``window_matcher.PathPredicate`` would replace a guard
  the window does not hold; ``windows.py`` writes no bound for any such
  guard, so a participant whose window is disconnected (the only case a path
  guard would be needed) is refused by name (``refused:path_predicate``), and
  a context template whose port reaches a physical end outside its window
  (the context differs from the base by a read beyond the window: a
  whole-species guard) is refused (``refused:whole_species_guard``).

Count and weight (item 2)
-------------------------
On a witness, ``T_ch`` is the number of matcher transitions of the
template(s) whose mapped product key (``rmg_moves.product_key`` of
``window_matcher.product_rows``) is a key of the channel; independently
``rmg_moves.generate_moves`` on the same molecules (participants as separate
components; an A + A record as two copies) counts the directed moves of the
channel (allowed status, the record's template, the target's product
species) and converts them by the slot-order rule (``SLOT_ORDER_RULE``).
Unequal counts refuse the template with both numbers
(``refused:count_mismatch``).  On the record's own witness (its reactant
graphs) the record's base-context templates jointly give ``T_R`` and the
per-transition weight ``w = ssa_multiplier * k / T_R`` over the whole channel
with the certified production rate (``k_table``) and its provenance
(estimate or library).  ``T_R = 0`` or a non-finite or non-positive ``w``
raises ``EmissionError``.  The calibration identity ``sum(w over the T_R
transitions) = ssa * k`` is checked at ``CALIBRATION_TEMPERATURES``
(necessary, not sufficient).

The record's graphs are a proxy whose artificial ends stand for the chain
beyond.  A base template whose window or port exteriors reach the padding
the windows added there (``reaches_beyond``: atoms the unextended witness
lacks) cannot apply to the proxy; the record row and the unextended-context
row are then ``not_applicable`` (never a zero count), and each template's
weight and calibration come from its emitted witness (``weight_source:
emitted_witness``) through the same production path.

A transfer whose products are its reactants (an unchanged_species move:
a radical moves between equal chains) is dropped by RMG's generator
(``_create_reaction``) but still has a rate: ``ProductionRates`` then builds
the witness's own channel as one explicit labelled ``TemplateReaction``
(``explicit_reaction``: rmg_moves mappings, RMG's template labels, RMG's
degeneracy rule) and sends it through the same ``_linked_family_pair``, so
the family's estimator gives k (provenance ``path: explicit_identity``).

Locality (item 3)
-----------------
Witnesses (``locality_witnesses``): the emitted one (the context's
fixed-point witness) plus extension vectors from ``boundaries`` (every
artificial port alone at +1 and +2, all ports jointly at +1 and +2, the
ports of each participant at +1 and +2), deduplicated as marked graphs
(``boundaries.marked_duplicates``; distinctness of the kept ones is
asserted), plus the context's own unextended witness.  On each the target is
that witness's own mapped product; ``k`` comes from the production path
(``ProductionRates``: the compiler's own generation and
``_linked_family_pair``; the production selection is rebound to the
record's: source, entry, rank, template direction and reference must match,
never re-estimated, no fallback); ``T_ch`` is recounted with BOTH the matcher
and ``rmg_moves``.  The record's own weight is a row too.  A template is
accepted iff every row is ``ok`` and ``max - min`` of ``log10 w`` over all rows
is at most ``LOCALITY_TOLERANCE_LOG10`` at every grid temperature; else
``refused:non_local_kinetics`` with the table.

Ladder (item 4)
---------------
``whole_species_cut`` cuts a ladder member's whole species (every atom
closed); no locality test.  The DP L + 1 control is checked by
``ladder_control``: its member's template is local (has a port) and covers
its own event.

Policy for RMG-filtered moves (item 5)
--------------------------------------
``compare_reference``: ``forbidden`` moves are excluded (a matcher key RMG
files only as forbidden fails); ``unchanged_species`` moves are required
(missing fails) and reported as matcher-only by policy.

Catalog collisions
------------------
Templates of different records with equal executable keys would each fire
on every matched transition.  ``resolve_collisions`` (one catalog) merges
copies whose weights agree to ``COLLISION_RTOL`` at every grid temperature
into one template listing every source record (the others ``merged``) and
refuses unequal ones (``refused:collision``: keys, records, weights).
``loaded_templates`` is what a catalog loads; the reference comparison fails
when a transition is produced more than once.

Schema (item 6)
---------------
``executable_payload`` / ``validate_payload``; ``attach`` sets
``executable_template`` on an ``EventRecord`` (schema ``rewrite_record/0.2``;
the template is not part of the event id, it carries its own sha256; only
the record's own payload attaches: ``ownership_problem``);
``matcher_from_payload`` is the temperature-specific projection to
``window_matcher.Template``.

Determinism (item 7)
--------------------
Canonical atom order and participant numbering; JSON with sorted keys;
witness atom ids appear only in the ``maps`` section, so the executable
section is byte-identical across processes, hash seeds and permuted input
atom orders.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable, Optional, Sequence

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc import window_matcher as wm

SCHEMA = "kmc-executable-template/1"
RECORD_SCHEMA = "rewrite_record/0.2"
LEGACY_RECORD_SCHEMA = "rewrite_record/0.1"
LOCALITY_TOLERANCE_LOG10 = 0.01
CALIBRATION_TEMPERATURES = (600.0, 700.0, 800.0)
CALIBRATION_RTOL = 1e-9
EXTENSION_STEPS = (1, 2)
CANONICAL_LEAF_LIMIT = 20000
BOUNDARY_LABELS = ("artificial", "physical", "interior")
TEMPLATE_STATUSES = ("emitted", "merged", "refused")
RECORD_STATUSES = ("emitted", "partial", "refused")
ESTIMATE_KINDS = ("RMG family estimate", "reference-thermo reverse",
                  "synthetic:unchanged_species")
SLOT_ORDER_RULE = (
    "one transition per distinct mapped outcome (rmg_moves.product_key): every "
    "directed mapping reaching one mapped product -- both template slot orders of "
    "a pair when the template is symmetric (identical or unequal radicals "
    "recombining: one product per site pair), every reactant resonance form, every "
    "subgraph isomorphism -- is one transition; an ordered template "
    "(disproportionation: acceptor/donor) keeps the two slot orders apart because "
    "they give different mapped products; an A + A pair is two separate copies")
NAMESPACES = {
    "executable": "template atom index: window_matcher.Template atom order (canonical)",
    "template": "RMG family template labels of the representative mapping",
    "witness": "emission witness atom ids (boundaries.Witness, explicit H)",
    "record": "the record's reactant-graph atom positions (certificates.witness_molecules)",
}


class EmissionError(ValueError):
    """A compile error: emission cannot proceed (``T_ch = 0``, a non-finite or
    non-positive weight, an edit on a port, a malformed payload)."""


def canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False)


def sha256_json(value) -> str:
    return hashlib.sha256(canonical_json(value).encode("ascii")).hexdigest()


# ---------------------------------------------------------------------------
# rows, states, molecules
# ---------------------------------------------------------------------------
def template_sides(witness, reverse: bool = False):
    """(reactant rows, product rows) of the template direction (dicts)."""
    r, p = dict(witness.reactant_rows()), dict(witness.product_rows())
    return (p, r) if reverse else (r, p)


def melt_state(rows):
    """A MeltState holding ``rows`` (atom id -> rmg_moves row), one component
    per connected part, atom ids kept."""
    from rmgpy.kmc.melt_state import MeltState, Node, Origin
    nodes, origin = {}, {}
    for k, part in enumerate(rm.split_rows(dict(rows))):
        for a, row in part.items():
            nodes[a] = Node(row[0], row[1], int(row[2]), int(row[3]), int(row[4]),
                            {n: float(o) for n, o in row[5]})
            origin[a] = Origin("chain", k, "witness", str(a), False)
    return MeltState(nodes, origin)


def first_form_rows(rows) -> dict:
    """``rows`` with every component written in the matcher's first resonance
    form (``window_matcher.component_forms``)."""
    state = melt_state(rows)
    out = {}
    for cid in state.components.ids():
        out.update(wm.form_rows(state, cid, wm.component_forms(state, cid)[0]))
    return out


def molecules_of(rows) -> list:
    """RMG molecules (atom ids kept), one per component, by least atom id."""
    return [rm.molecule_from_rows(part) for part in rm.split_rows(dict(rows))]


def _cyclic(rows) -> set:
    from rmgpy.kmc import windows as win
    return win._cyclic_atoms(rows)


def _witness_ports(witness) -> dict:
    """atom id -> (participant, port id, tag, physicality) for every witness
    port atom and declared cap (construction metadata)."""
    out = {}
    for participant in witness.participants:
        for port in participant.metadata.ports:
            info = (participant.index, port.port_id, port.tag, port.physicality)
            out[port.atom_id] = info
            for cap in port.cap_atom_ids:
                out[cap] = info
    return out


# ---------------------------------------------------------------------------
# canonical labelling
# ---------------------------------------------------------------------------
def _canonical_order(n, exec_labels, meta_labels, edges):
    """The atom order minimizing (executable serialization, metadata
    serialization) over individualization-refinement (``edges``:
    {(i, j): label} with string labels).  Interchangeable hydrogen leaves
    (degree one, same neighbour and bond, same labels) are automorphic and
    are not branched."""
    nbrs = [[] for _ in range(n)]
    for (i, j), lab in edges.items():
        nbrs[i].append((j, lab))
        nbrs[j].append((i, lab))

    def rank(keys):
        order = {k: r for r, k in enumerate(sorted(set(keys)))}
        return [order[k] for k in keys]

    def refine(colors):
        count = len(set(colors))
        while True:
            new = rank([(colors[v], tuple(sorted((lab, colors[u]) for u, lab in nbrs[v])))
                        for v in range(n)])
            c = len(set(new))
            if c == count:
                return new
            colors, count = new, c

    def serialize(order):
        pos = {v: k for k, v in enumerate(order)}
        ex = tuple(exec_labels[v] for v in order)
        ed = tuple(sorted((min(pos[i], pos[j]), max(pos[i], pos[j]), lab)
                          for (i, j), lab in edges.items()))
        me = tuple(meta_labels[v] for v in order)
        return (ex, ed), me

    best = []
    leaves = [0]

    def search(colors):
        cells = defaultdict(list)
        for v, c in enumerate(colors):
            cells[c].append(v)
        target = next((cells[c] for c in sorted(cells) if len(cells[c]) > 1), None)
        if target is None:
            leaves[0] += 1
            if leaves[0] > CANONICAL_LEAF_LIMIT:
                raise EmissionError(f"canonical labelling exceeds {CANONICAL_LEAF_LIMIT} leaves")
            order = sorted(range(n), key=lambda v: colors[v])
            ser = serialize(order)
            if not best or ser < best[0][0]:
                best[:] = [(ser, order)]
            return
        if all(len(nbrs[v]) == 1 for v in target) and len({nbrs[v][0] for v in target}) == 1:
            choices = [min(target)]   # automorphic leaves: one branch serializes them all
        else:
            choices = sorted(target)
        for v in choices:
            search(refine(rank([(c, 0 if u == v else 1) for u, c in enumerate(colors)])))

    search(refine(rank([(exec_labels[v], meta_labels[v]) for v in range(n)])))
    return best[0][1], best[0][0]


# ---------------------------------------------------------------------------
# item 1: the template cut from a window
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CutTemplate:
    """One executable template in canonical atom order.

    ``atoms``: per atom (participant, element, isotope, charge, radical,
    lone pairs, kind, exterior, final attributes); ``bonds``: (i, j, order
    before, order after) with i < j (0: absent on that side); ``meta``: per
    atom (boundary label or None, RMG labels, witness port or None);
    ``witness_atoms``: the emission witness atom id per atom (the full-H map);
    ``components``: per participant the least witness atom id of its
    template-direction reactant component; ``key``: sha256 of the executable
    content (atoms, bonds)."""

    family: str
    direction: str
    pathway: tuple
    atoms: tuple
    bonds: tuple
    meta: tuple
    witness_atoms: tuple
    components: tuple
    context: tuple = ()
    key: str = ""

    @property
    def participants(self) -> int:
        return len(self.components)

    @property
    def ports(self) -> tuple:
        return tuple(i for i, a in enumerate(self.atoms) if a[6] == "port")

    @property
    def whole_species(self) -> bool:
        return not self.ports

    def boundary_labels(self) -> dict:
        return {i: self.meta[i][0] for i in self.ports}

    def edited(self) -> tuple:
        out = {i for i, a in enumerate(self.atoms) if a[8]}
        out |= {x for i, j, before, after in self.bonds if before != after for x in (i, j)}
        return tuple(sorted(out))

    def participant_connected(self) -> bool:
        """Whether every participant's window is one connected piece."""
        adj = defaultdict(set)
        for i, j, before, _after in self.bonds:
            if before:
                adj[i].add(j)
                adj[j].add(i)
        for p in range(self.participants):
            members = [i for i, a in enumerate(self.atoms) if a[0] == p]
            seen, stack = {members[0]}, [members[0]]
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            if set(members) - seen:
                return False
        return True

    def executable(self) -> dict:
        """The id-free executable content (what the hash covers)."""
        return {"atoms": [_atom_json(a) for a in self.atoms],
                "bonds": [[i, j, before, after] for i, j, before, after in self.bonds],
                "participants": self.participants}

    def matcher(self, weight: float, name: Optional[str] = None) -> wm.Template:
        """The ``window_matcher.Template`` of this cut at one weight."""
        atoms = [wm.WindowAtom(element=a[1], participant=a[0], kind=a[6], exterior=a[7],
                               isotope=a[2], charge=a[3], radical=a[4], lone_pairs=a[5])
                 for a in self.atoms]
        bonds = [(i, j, before) for i, j, before, _after in self.bonds if before]
        edits = [wm.BondEdit(i, j, after) for i, j, before, after in self.bonds
                 if before != after]
        edits += [wm.AtomEdit(i, **dict(a[8])) for i, a in enumerate(self.atoms) if a[8]]
        return wm.Template(name or self.name, self.family, atoms, bonds, edits, weight=weight)

    @property
    def name(self) -> str:
        return f"{self.family}:{self.direction}:{self.key[:16]}"


def _atom_json(a) -> dict:
    return {"participant": a[0], "element": a[1], "isotope": a[2], "charge": a[3],
            "radical": a[4], "lone_pairs": a[5], "kind": a[6],
            "exterior": [list(x) for x in a[7]], "final": {k: v for k, v in a[8]}}


def cut_template(witness, atoms, *, family: str, pathway, direction: str,
                 reverse: bool = False, context=(), labels=None) -> CutTemplate:
    """Cut the window ``atoms`` of ``witness`` as an executable template (see
    the module docstring, item 1).  ``labels``: witness atom id -> RMG labels
    (template namespace metadata)."""
    window = set(atoms)
    r0, p0 = template_sides(witness, reverse)
    unknown = sorted(window - set(r0))
    if unknown:
        raise EmissionError(f"window atoms {unknown[:5]} are not witness atoms")
    try:
        r, p = first_form_rows(r0), first_form_rows(p0)
    except wm.WindowMatcherError as error:
        raise EmissionError(f"resonance normalisation of the witness failed: {error}") from None
    comps = rm.split_rows(r)
    comp_of = {a: k for k, part in enumerate(comps) for a in part}
    used = sorted({comp_of[a] for a in window})
    if len(used) not in (1, 2):
        raise EmissionError(f"the window touches {len(used)} reactant components (1 or 2)")
    cyclic = _cyclic(r)
    ports_of = _witness_ports(witness)
    labels = labels or {}
    order_atoms = sorted(window)
    raw = {}
    for a in order_atoms:
        out_r = sorted((n, float(o)) for n, o in r[a][5] if n not in window)
        out_p = sorted((n, float(o)) for n, o in p[a][5] if n not in window)
        if out_r != out_p:
            raise EmissionError(f"an edit crosses the window boundary at atom {a} "
                                f"({out_r} -> {out_p})")
        kind = "port" if out_r else "closed"
        exterior = tuple(sorted((o, r[n][0], n in cyclic) for n, o in out_r))
        final = tuple((name, p[a][k]) for name, k in (("charge", 2), ("radical", 3),
                                                      ("lone_pairs", 4)) if r[a][k] != p[a][k])
        raw[a] = (used.index(comp_of[a]), r[a][0], r[a][1], int(r[a][2]), int(r[a][3]),
                  int(r[a][4]), kind, exterior, final)
    bonds = {}
    for side, rows in ((0, r), (1, p)):
        for a in order_atoms:
            for n, o in rows[a][5]:
                if n in window and a < n:
                    bonds.setdefault((a, n), [0.0, 0.0])[side] = float(o)
    edited = {a for a in order_atoms if raw[a][8]} | {
        x for pair, (before, after) in bonds.items() if before != after for x in pair}
    on_port = sorted(a for a in edited if raw[a][6] == "port")
    if on_port:
        raise EmissionError(f"an edit touches the window port(s) {on_port}: edits never "
                            "touch a port")
    # boundary labels from the construction metadata: the witness ports the
    # exterior of each port atom reaches outside the window
    reach_cache = {}

    def reach(start):
        if start in reach_cache:
            return reach_cache[start]
        seen, stack, found = {start}, [start], set()
        while stack:
            x = stack.pop()
            if x in ports_of:
                found.add(ports_of[x])
            for n, _o in r[x][5]:
                if n not in window and n not in seen:
                    seen.add(n)
                    stack.append(n)
        value = frozenset(found)
        for x in seen:
            reach_cache[x] = value
        return value

    meta = {}
    for a in order_atoms:
        boundary = None
        if raw[a][6] == "port":
            found = set()
            for n, _o in r[a][5]:
                if n not in window:
                    found |= reach(n)
            phys = {x[3] for x in found}
            boundary = "artificial" if "artificial" in phys else \
                "physical" if "physical" in phys else "interior"
        wp = ports_of.get(a)
        meta[a] = (boundary, tuple(sorted(labels.get(a, ()))),
                   None if wp is None else (wp[1], wp[2], wp[3]))
    perms = [(0,), ] if len(used) == 1 else [(0, 1), (1, 0)]
    edge_labels = {pair: json.dumps(v) for pair, v in bonds.items()}
    best = None
    index = {a: k for k, a in enumerate(order_atoms)}
    n = len(order_atoms)
    edges = {(index[x], index[y]): lab for (x, y), lab in edge_labels.items()}
    for perm in perms:
        exec_labels = [json.dumps(_atom_json((perm[raw[a][0]],) + raw[a][1:]), sort_keys=True)
                       for a in order_atoms]
        meta_labels = [json.dumps(meta[a]) for a in order_atoms]
        order, ser = _canonical_order(n, exec_labels, meta_labels, edges)
        if best is None or ser < best[0]:
            best = (ser, perm, order)
    _ser, perm, order = best
    ordered = [order_atoms[k] for k in order]
    pos = {a: k for k, a in enumerate(ordered)}
    atoms_out = tuple((perm[raw[a][0]],) + raw[a][1:] for a in ordered)
    bonds_out = tuple(sorted((min(pos[x], pos[y]), max(pos[x], pos[y]), v[0], v[1])
                             for (x, y), v in bonds.items()))
    components = [None] * len(used)
    for k, c in enumerate(used):
        components[perm[k]] = min(comps[c])
    cut = CutTemplate(family, direction, tuple(pathway or ()), atoms_out, bonds_out,
                      tuple(meta[a] for a in ordered), tuple(ordered), tuple(components),
                      tuple(context))
    key = sha256_json(cut.executable())
    return CutTemplate(family, direction, cut.pathway, atoms_out, bonds_out, cut.meta,
                       cut.witness_atoms, cut.components, cut.context, key)


def whole_species_cut(witness, *, family: str, pathway, direction: str, reverse=False,
                      labels=None) -> CutTemplate:
    """Item 4: a ladder member's whole species as one template (every atom of
    the reacting components closed; matches that species only)."""
    r, _p = template_sides(witness, reverse)
    return cut_template(witness, set(r), family=family, pathway=pathway, direction=direction,
                        reverse=reverse, labels=labels)


def labels_of(move, reverse: bool = False) -> dict:
    """witness atom id -> RMG labels of a ``rmg_moves`` move on the template
    direction's reactant side (the move's reactant labels, or its product
    labels when the template runs opposite to the move)."""
    out = defaultdict(set)
    items = ((lab, a) for lab, a in (move.product_labels or ())) if reverse else \
        ((lab, a) for _s, lab, a in move.labels)
    for lab, a in items:
        out[a].add(lab)
    return {a: tuple(sorted(v)) for a, v in out.items()}


def template_refusal(cut: CutTemplate) -> Optional[str]:
    """Named refusals of a cut (module docstring, item 1)."""
    if not cut.participant_connected():
        return ("refused:path_predicate: a participant's window is disconnected; a path "
                "guard would replace a read windows.py gives no written bound for")
    physical = sorted(i for i, label in cut.boundary_labels().items() if label == "physical")
    if physical:
        return (f"refused:whole_species_guard: port atom(s) {physical} reach only a physical "
                "end outside the window (the context differs by a read beyond the window)")
    return None


# ---------------------------------------------------------------------------
# item 2: counts on one witness
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RmgChannel:
    """RMG's channel on one witness: the distinct mapped outcomes (``keys``)
    of the allowed directed moves with the channel's template and product
    species; ``directed`` raw mappings; ``slot_orders`` per outcome."""

    keys: tuple
    directed: int
    slot_orders: tuple
    statuses: tuple

    @property
    def transitions(self) -> int:
        return len(self.keys)


class RmgCounter:
    """The independent count: ``rmg_moves.generate_moves`` (cached per
    witness signature) and species keys."""

    def __init__(self, families, resonance: str = "all", max_entries: int = 64):
        from rmgpy.kmc import channels as chm
        self.families = families
        self.resonance = resonance
        self.species = chm.SpeciesKeys()
        self._moves = {}
        self._max = max_entries
        self._product_species = {}

    def family(self, name):
        return self.families[name] if isinstance(name, str) else name

    def moves(self, family, molecules, direction):
        family = self.family(family)
        if direction == "reverse" and family.reverse_template is None:
            direction = "forward"
        key = (family.label, direction, tuple(rm._signature(rm.molecule_rows(m))
                                              for m in molecules))
        if key not in self._moves:
            if len(self._moves) >= self._max:
                self._moves.pop(next(iter(self._moves)))
            self._moves[key] = rm.generate_moves(family, list(molecules), direction,
                                                 resonance=self.resonance)
        return self._moves[key]

    def species_of_rows(self, rows) -> tuple:
        sig = rm._signature(dict(rows))
        if sig not in self._product_species:
            self._product_species[sig] = tuple(sorted(
                self.species(rm.molecule_from_rows(part)) for part in rm.split_rows(dict(rows))))
        return self._product_species[sig]

    def channel(self, family, molecules, direction, pathway, product_species) -> RmgChannel:
        keys, orders, statuses = {}, defaultdict(set), Counter()
        directed = 0
        for m in self.moves(family, molecules, direction):
            statuses[m.status] += 1
            if m.status not in ("kept", "unchanged_species") or m.products is None:
                continue
            if m.template is None:
                continue
            if pathway is not None and tuple(m.template) != tuple(pathway):
                continue
            if self.species_of_rows(m.product_rows()) != tuple(product_species):
                continue
            directed += 1
            k = rm.product_key(m)
            keys[k] = True
            orders[k].add(m.order)
        ordered = tuple(sorted(keys))
        return RmgChannel(ordered, directed,
                          tuple(tuple(sorted(orders[k])) for k in ordered),
                          tuple(sorted(statuses.items())))


@dataclass(frozen=True)
class Count:
    matcher: int
    rmg: int
    rmg_directed: int
    target_found: bool
    foreign: int
    rule: str = SLOT_ORDER_RULE

    @property
    def equal(self) -> bool:
        return self.matcher == self.rmg

    def as_dict(self) -> dict:
        return {"matcher": self.matcher, "rmg": self.rmg, "rmg_directed": self.rmg_directed,
                "target_found": self.target_found, "foreign": self.foreign}


def count_pathway(family, pathway, reverse: bool):
    """The pathway filter of the channel count: the record's template, except
    for a family that is its own reverse read opposite to its template side
    (RMG reads another node on those reactants: product species decide)."""
    fam = family
    if getattr(fam, "reverse_template", "x") is None and reverse:
        return None
    return tuple(pathway) if pathway else None


def count_on(templates: Sequence[wm.Template], reactant_rows, product_rows, *, family,
             direction, pathway, counter, product_species=None, target=True) -> tuple:
    """(Count, transitions in the channel, channel) of ``templates`` on one
    witness (template-direction rows; components separate).  ``target``: the
    product rows are the witness's own mapped product (its key must be one of
    the transitions)."""
    species = product_species if product_species is not None else \
        counter.species_of_rows(product_rows)
    channel = counter.channel(family, molecules_of(reactant_rows), direction, pathway, species)
    state = melt_state(reactant_rows)
    result = wm.match(list(templates), state)
    keyset = set(channel.keys)
    inside, foreign = [], 0
    for tr in result.transitions:
        key = rm.product_key(wm.product_rows(state, tr))
        if key in keyset:
            inside.append((tr, key))
        else:
            foreign += 1
    found = True
    if target:
        found = rm.product_key(dict(product_rows)) in {k for _t, k in inside}
    return (Count(len(inside), channel.transitions, channel.directed, found, foreign),
            inside, channel)


# ---------------------------------------------------------------------------
# rates and weights
# ---------------------------------------------------------------------------
def rate_kind(rate_source: dict) -> str:
    """``estimate`` (an RMG selection, or its reference-thermo reverse) or
    ``library`` (anything else: library, archived or attached rates)."""
    kind = (rate_source or {}).get("kind")
    return "estimate" if kind in ESTIMATE_KINDS else "library"


def weight_table(temperatures, k, ssa_multiplier, transitions) -> list:
    """``w = ssa * k / T_ch`` per grid temperature (``EmissionError`` on
    ``T_ch = 0`` or a non-finite or non-positive weight)."""
    if not transitions or transitions <= 0:
        raise EmissionError("T_ch = 0: the template makes no transition of its channel")
    out = []
    for t, value in zip(temperatures, k):
        w = float(ssa_multiplier) * float(value) / int(transitions)
        if not math.isfinite(w) or w <= 0:
            raise EmissionError(f"non-finite or non-positive weight {w!r} at {t} K")
        out.append(w)
    return out


def calibration(temperatures, k, ssa_multiplier, weights, transitions) -> list:
    """The identity ``sum(w over T_ch transitions) = ssa * k`` at the
    calibration temperatures (rows: T, sum, ssa*k, relative error, ok)."""
    rows = []
    at = {float(t): i for i, t in enumerate(temperatures)}
    for t in CALIBRATION_TEMPERATURES:
        i = at.get(float(t))
        if i is None:
            rows.append({"T": t, "ok": False, "reason": "not a grid temperature"})
            continue
        total = math.fsum([weights[i]] * int(transitions))
        want = float(ssa_multiplier) * float(k[i])
        err = abs(total - want) / want
        rows.append({"T": t, "sum_w": total, "ssa_k": want, "relative_error": err,
                     "ok": err <= CALIBRATION_RTOL})
    return rows


@dataclass(frozen=True)
class RateRequest:
    """One production-rate question: the channel of ``witness``'s own mapped
    event, generated on the case-direction reactants (``generation``) with
    the case-direction products (``products``); the record's ``template``
    and ``orientation``; the original ``rate_source`` to rebind to."""

    family: str
    template: str
    orientation: str
    generation: Any
    products: Any
    witness: Any = None
    label: str = ""
    vector: tuple = ()
    rate_source: Any = None
    owner_rate_source: Any = None


@dataclass(frozen=True)
class RateResult:
    temperatures: tuple
    k: tuple
    ssa_multiplier: float
    provenance: dict = field(default_factory=dict)
    problems: tuple = ()


_SELECTION_KEYS = ("kind", "reference", "source", "entry", "rank", "family_template_direction")


def _selection(rate_source) -> dict:
    rs = rate_source or {}
    out = {k: rs.get(k) for k in _SELECTION_KEYS if k in rs}
    forward = rs.get("forward_rate_source")
    if isinstance(forward, dict):
        out["forward"] = {k: forward.get(k) for k in _SELECTION_KEYS if k in forward}
    return out


def rebind_problems(original, new) -> list:
    """The production selection on a witness against the record's: every
    selection field (source, entry, rank, template direction, reference) must
    match (never re-estimated)."""
    a, c = _selection(original), _selection(new)
    out = []
    for key in sorted(set(a) | set(c)):
        if a.get(key) != c.get(key):
            out.append(f"rate source {key}: record {a.get(key)!r} vs witness {c.get(key)!r}")
    return out


class ProductionRates:
    """``k`` on a witness through the compiler's own production path: RMG's
    generation (``generate_reactions_from_families`` on the case-direction
    reactants, restricted to the witness's products), then the compiler's
    ``_linked_family_pair`` (estimate, thermo conversion, reference-thermo
    reverse), then the record with the requested orientation and template.
    The selection is rebound to the record's (``rebind_problems``) and, for
    a reference-thermo reverse, the pair's estimating record to the record's
    rate owner.  Cached per (family, generation signature, products,
    template)."""

    def __init__(self, database, families, *, database_path=None, database_sha=None,
                 temperature_grid=None, max_entries: int = 256):
        from rmgpy.kmc.compiler import DEFAULT_T_GRID, EventSetCompiler
        self.database = database
        self.compiler = EventSetCompiler(
            database.kinetics, (), tuple(families), database_path=database_path,
            thermo_database=database.thermo, rmg_database_sha=database_sha,
            temperature_grid=tuple(temperature_grid or DEFAULT_T_GRID))
        self._cache = {}
        self._max = max_entries
        self.calls = 0
        self.generated = 0
        self.explicit_calls = 0
        self._counter = None

    def counter(self):
        if self._counter is None:
            self._counter = RmgCounter(self.database.kinetics.families)
        return self._counter

    def _linked(self, species, reaction):
        from rmgpy.kmc.compiler import SiteProxy
        proxy = SiteProxy("i138_emission_witness", tuple(species),
                          participant_site_types=tuple(f"w{k}" for k in range(len(species))),
                          metadata={"centred": True, "context_class": "emission_witness",
                                    "proxy_units": 0, "family_candidates": [reaction.family]})
        try:
            records, _initiating = self.compiler._linked_family_pair(proxy, reaction, {})
            return tuple(records), ""
        except Exception as error:   # reported, never hidden
            return None, f"production pair failed: {type(error).__name__}: {error}"

    def explicit_reaction(self, request):
        """(reaction, problem): the witness's own channel as one explicit
        labelled ``TemplateReaction`` in the family's template direction
        (``request.generation`` is that side), built without RMG's generator,
        which drops a reaction whose products are its reactants.

        The channel's moves are ``rmg_moves`` mappings (RMG's own template
        matching and recipe; template labels by ``get_reaction_template_labels``,
        for an identity transfer through ``_unchanged_species_template``).  The
        degeneracy is RMG's: the distinct mapped outcomes of the channel (the
        generator merges isomorphic, non-identical reactions and drops identical
        ones, ``find_degenerate_reactions``), halved for two identical reactants
        (``reduce_same_reactant_degeneracy``).  Kinetics then come from the
        family's estimator inside the compiler's ``_linked_family_pair``, the
        same call the generator's reaction goes through."""
        from rmgpy.data.kinetics.family import TemplateReaction
        from rmgpy.species import Species
        fam = self.database.kinetics.families[request.family]
        counter = self.counter()
        gen = molecules_of(request.generation)
        species = counter.species_of_rows(request.products)
        pathway = tuple(request.template.split(";")) if request.template else None
        channel = counter.channel(fam, gen, "forward", pathway, species)
        if not channel.transitions:
            return None, ("explicit path: no forward rmg_moves of the template reach the "
                          "witness's products (the generation side is not the family's "
                          "template side)")
        degeneracy = float(channel.transitions)
        if len(gen) == 2 and gen[0].is_isomorphic(gen[1]):
            degeneracy *= 0.5
        reactants = [Species(molecule=[m.copy(deep=True)]) for m in gen]
        products = [Species(molecule=[m.copy(deep=True)]) for m in molecules_of(request.products)]
        for spc in reactants + products:
            spc.generate_resonance_structures(keep_isomorphic=False)
        reaction = TemplateReaction(reactants=reactants, products=products,
                                    degeneracy=degeneracy, reversible=fam.reversible,
                                    family=fam.label, is_forward=True,
                                    template=list(pathway or ()), electrons=fam.electrons)
        return (reaction, species, reactants), ""

    def explicit(self, request):
        """Records of the explicit labelled-reaction path (``explicit_reaction``)."""
        self.explicit_calls += 1
        built, why = self.explicit_reaction(request)
        if built is None:
            return None, why
        reaction, _species, reactants = built
        return self._linked(reactants, reaction)

    def _pair(self, request):
        from rmgpy.kmc.compiler import SiteProxy, _resolve_generator, _template
        from rmgpy.species import Species
        gen = molecules_of(request.generation)
        prods = molecules_of(request.products)
        key = (request.family, tuple(rm._signature(rm.molecule_rows(m)) for m in gen),
               tuple(sorted(rm._signature(rm.molecule_rows(m)) for m in prods)),
               request.template)
        if key in self._cache:
            return self._cache[key]
        self.generated += 1
        species = [Species(molecule=[m]) for m in gen]
        products = [Species(molecule=[m]) for m in prods]
        generator = _resolve_generator(self.compiler.kinetics_database,
                                       self.compiler.generation_provider)
        reactions = [x for x in generator([s.copy(deep=True) for s in species],
                                          products=[s.copy(deep=True) for s in products],
                                          only_families=[request.family], resonance=True)
                     if _template(x) == request.template]
        if not reactions and self.identity(request):
            # RMG's generator drops a transfer whose products are its reactants
            # (_create_reaction); it still has a rate: the explicit path
            records, why = self.explicit(request)
            value = (records, why, "explicit_identity")
        elif len(reactions) != 1:
            value = (None, f"production generation found {len(reactions)} reactions of the "
                           f"template on the witness (one needed)", "generator")
        else:
            records, why = self._linked(species, reactions[0])
            value = (records, why, "generator")
        if len(self._cache) >= self._max:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = value
        return value

    def identity(self, request) -> bool:
        """Whether the witness's products are its reactants (as species)."""
        counter = self.counter()
        return counter.species_of_rows(request.generation) == \
            counter.species_of_rows(request.products)

    def __call__(self, request: RateRequest) -> RateResult:
        self.calls += 1
        records, why, path = self._pair(request)
        if records is None:
            return RateResult((), (), 0.0, {}, (why,))
        hits = [x for x in records if x.orientation == request.orientation and
                x.template == request.template and x.k_table]
        if len(hits) != 1:
            return RateResult((), (), 0.0, {}, (f"production pair has {len(hits)} records of "
                                                 f"orientation {request.orientation}",))
        rec = hits[0]
        problems = list(rebind_problems(request.rate_source, rec.rate_source)) \
            if request.rate_source is not None else []
        if request.owner_rate_source is not None:
            owners = [x for x in records if x is not rec and x.template == request.template
                      and (x.rate_source or {}).get("kind") == "RMG family estimate"]
            if len(owners) != 1:
                problems.append(f"production pair has {len(owners)} estimating records")
            else:
                problems += [f"rate owner: {p}" for p in
                             rebind_problems(request.owner_rate_source, owners[0].rate_source)]
        problems = tuple(problems)
        return RateResult(tuple(rec.k_table["T"]), tuple(rec.k_table["k"]),
                          float(rec.ssa_multiplier),
                          {"selection": _selection(rec.rate_source), "event_id": rec.event_id,
                           "status": rec.status, "path": path}, problems)


# ---------------------------------------------------------------------------
# item 3: locality over extension vectors
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LocalityWitness:
    label: str
    vector: tuple
    status: str
    reason: str = ""
    witness: Any = field(default=None, compare=False)


def locality_plan(witness) -> list:
    """(label, vector) of the extension vectors: every artificial port alone at
    each step of ``EXTENSION_STEPS``, all ports jointly, the ports of each
    participant; the zero vector (the emitted witness) first."""
    ports = list(b.vector_ports(witness))
    n = len(ports)
    plan = [("emitted", (0,) * n)]
    for k, p in enumerate(ports):
        for s in EXTENSION_STEPS:
            v = [0] * n
            v[k] = s
            plan.append((f"port {p[0]}:{p[1]} +{s}", tuple(v)))
    if n:
        for s in EXTENSION_STEPS:
            plan.append((f"joint +{s}", (s,) * n))
    for q in sorted({p[0] for p in ports}):
        for s in EXTENSION_STEPS:
            plan.append((f"participant {q} +{s}",
                         tuple(s if p[0] == q else 0 for p in ports)))
    return plan


def locality_witnesses(witness, providers=None, *, ceiling=None, plan=None) -> list:
    """The locality witnesses (``locality_plan`` built by
    ``boundaries.evaluate_vector``), equal vectors folded, the ok ones
    deduplicated as marked graphs (``boundaries.marked_duplicates``) with
    the kept ones asserted distinct."""
    out, seen = [], {}
    for label, vector in (plan if plan is not None else locality_plan(witness)):
        if vector in seen:
            continue
        seen[vector] = label
        if not any(vector):
            out.append(LocalityWitness(label, vector, "ok", "", witness))
            continue
        res = b.evaluate_vector(witness, vector, providers, ceiling=ceiling)
        out.append(LocalityWitness(label, vector, res.status, res.reason, res.witness))
    ok = [(x.label, x.witness) for x in out if x.status == "ok"]
    dropped = {second for _first, second in b.marked_duplicates(ok)}
    kept = [x for x in out if x.label not in dropped]
    if b.marked_duplicates([(x.label, x.witness) for x in kept if x.status == "ok"]):
        raise EmissionError("locality witnesses are not distinct after deduplication")
    return kept


def spread(rows, temperatures) -> tuple:
    """(per-temperature max - min of log10 w over the rows, max)."""
    per = []
    for i, _t in enumerate(temperatures):
        values = [math.log10(r["w"][i]) for r in rows]
        per.append(max(values) - min(values) if values else float("inf"))
    return per, (max(per) if per else float("inf"))


# ---------------------------------------------------------------------------
# one candidate template of a record
# ---------------------------------------------------------------------------
def status_of(text: str) -> str:
    """``refused:<name>`` (or ``error:<name>``) of a refusal text."""
    parts = (text or "error").split(":")
    return ":".join(x.strip() for x in parts[:2]) if len(parts) > 1 else parts[0]


@dataclass(eq=False)
class Candidate:
    """A context representative's cut, with its sources (member label,
    context index and assignment, witness, window)."""

    cut: Optional[CutTemplate]
    member: str
    context_index: int
    assignment: tuple
    witness: Any
    window_atoms: frozenset
    certificate: Optional[str] = None
    base: bool = False
    refusal: Optional[str] = None
    unextended: Any = None
    sources: list = field(default_factory=list)


def context_candidates(case, results, merges, *, family, pathway, direction, reverse,
                       providers=None) -> tuple:
    """(candidates, context rows): one candidate per fixed-point context that
    merged into no earlier one (the base first), refused ones named; every
    other context's status reported."""
    from rmgpy.kmc import certificates as cm
    from rmgpy.kmc import windows as win
    merged = {j for _i, j in merges}
    into = {j: i for i, j in merges}
    labels = labels_of(case.moves[0], reverse)
    cands, rows = [], []
    for k, r in enumerate(results):
        assignment = tuple((tuple(p), e) for p, e in r.assignment)
        row = {"index": k, "assignment": [[list(p), e] for p, e in assignment if e],
               "status": r.status, "reason": (r.reason or "")[:240]}
        if r.status != "fixed_point":
            row["emission"] = "refused upstream (context status)"
            rows.append(row)
            continue
        if k in merged:
            row["emission"] = f"merged into context {into[k]}"
            rows.append(row)
            continue
        cert = None
        if r.observation is not None and r.observation.certificate is not None:
            cert = cm.certificate_sha256(r.observation.certificate)
        context = tuple((list(p), e) for p, e in assignment if e)
        try:
            cut = cut_template(r.witness, r.window.atoms, family=family, pathway=pathway,
                               direction=direction, reverse=reverse,
                               context=tuple((p[0], p[1], e) for p, e in assignment if e),
                               labels=labels)
            refusal = template_refusal(cut)
        except EmissionError as error:
            cut, refusal = None, f"error:{error}"
        unextended = case.witness
        if k:
            try:
                unextended = win.reconstruct(case.witness, assignment, providers)
            except b.BoundaryError as error:
                refusal = refusal or f"refused:end: {error}"
        cands.append(Candidate(cut, case.label, k, assignment, r.witness,
                               frozenset(r.window.atoms), cert, k == 0, refusal, unextended))
        row["emission"] = "candidate" if refusal is None else refusal
        row["context"] = context
        rows.append(row)
    return cands, rows


def reaches_beyond(cut, witness, unextended) -> bool:
    """Whether the cut's window atoms or its port exteriors (atoms of
    ``witness``) reach atoms the ``unextended`` witness lacks: the padding
    added at its artificial ends.  Such a template cannot apply to the
    unextended witness (nor to the record's proxy graph)."""
    if unextended is None or cut is None:
        return False
    r, p = dict(witness.reactant_rows()), dict(witness.product_rows())
    reach = set(cut.witness_atoms)
    for a in cut.ports:
        # a port is never edited (``cut_template``), so its neighbours are
        # the same on both sides; one side is read, the other checked
        aid = cut.witness_atoms[a]
        neighbours = {n for n, _o in r[aid][5]}
        if neighbours != {n for n, _o in p[aid][5]}:
            raise EmissionError(f"port atom {aid} has different neighbours on the two sides")
        reach |= neighbours
    have = set(dict(unextended.reactant_rows())) | set(dict(unextended.product_rows()))
    return not reach <= have


def dedupe(cands) -> list:
    """Candidates with equal executable keys are one template (sources kept)."""
    out, by = [], {}
    for c in cands:
        if c.cut is None:
            out.append(c)
            continue
        if c.cut.key in by:
            first = by[c.cut.key]
            first.sources.append((c.member, c.context_index))
            first.base = first.base or c.base
            continue
        c.sources.append((c.member, c.context_index))
        by[c.cut.key] = c
        out.append(c)
    return out


# ---------------------------------------------------------------------------
# records
# ---------------------------------------------------------------------------
def record_direction(record, family) -> str:
    """The family direction a record runs in (its template side is its
    reactants: forward; else reverse; an own-reverse family: forward)."""
    from rmgpy.kmc import certificates as cm
    if getattr(family, "reverse_template", "x") is None:
        return "forward"
    side, _a = cm._template_fields(_as_dict(record))
    return "forward" if side in ("reactant_graphs", "rate_witness_reactant_graphs") else "reverse"


def _as_dict(record) -> dict:
    return record.to_dict() if hasattr(record, "to_dict") else dict(record)


def record_graph_rows(record, which: str = "reactant_graphs") -> dict:
    """atom id -> row of the record's own witness (``witness_molecules`` ids,
    participants separate; an A + A record whose graphs list one species
    with multiplicity 2 is two copies)."""
    from rmgpy.kmc import certificates as cm
    record = _as_dict(record)
    graphs = list(record[which])
    if which == "reactant_graphs":
        mult = list(record.get("reactant_multiplicities") or [])
        if len(graphs) == 1 and mult == [2]:
            graphs = graphs * 2
    rows = {}
    for m in cm.witness_molecules(graphs):
        rows.update(rm.molecule_rows(cm.aromatic_root(m)))
    return rows


def upstream_refusal(record, record_certificate, by_id=None) -> Optional[str]:
    """The upstream status of a record (None: certified and firing)."""
    from rmgpy.kmc import certificates as cm
    record = _as_dict(record)
    if record.get("status") not in cm.FIRING_STATUSES:
        return f"refused:record_status: {record.get('status_reason', '')[:200]}"
    if record_certificate is None:
        return "refused:certificate: no record certificate"
    if by_id is not None and not cm.requires_certification(record, by_id):
        return ("refused:certificate: outside the required population "
                f"({(record_certificate.get('gaps') or [''])[0][:160]})")
    if not record_certificate.get("certified"):
        return f"refused:certificate: {(record_certificate.get('gaps') or [''])[0][:200]}"
    return None


def rate_of(record) -> dict:
    record = _as_dict(record)
    table = record.get("k_table") or {}
    rs = record.get("rate_source") or {}
    return {"T": list(table.get("T") or ()), "k": list(table.get("k") or ()),
            "ssa_multiplier": float(record.get("ssa_multiplier", 1.0)),
            "provenance": {"kind": rs.get("kind"), "class": rate_kind(rs),
                           "selection": _selection(rs), "event_id": record.get("event_id"),
                           "orientation": record.get("orientation")}}


def emit_record(record, *, family, members, counter, rates, record_certificate=None,
                by_id=None, providers=None, ceiling=None, check_upstream=True,
                owner_rate_source=None) -> dict:
    """The executable-template payload of one record.

    ``members``: [(case, run_member row, (context results, merges))] of the
    record's channel members on provider witnesses.  Upstream refusals
    (record status, certificate, window status) are carried through as
    refused; emission refusals are named per template."""
    record = _as_dict(record)
    fam = counter.family(family)
    direction = record_direction(record, fam)
    if owner_rate_source is None and by_id is not None and \
            (record.get("rate_source") or {}).get("kind") == "reference-thermo reverse":
        owner_rate_source = (by_id.get(record.get("reverse_of")) or {}).get("rate_source")
    pathway = tuple(record["template"].split(";")) if record.get("template") else ()
    rate = rate_of(record)
    head = {"schema": SCHEMA, "family": fam.label, "direction": direction,
            "pathway": list(pathway), "record_event_id": record.get("event_id"),
            "orientation": record.get("orientation"), "rate": rate,
            "certificate": (record.get("rate_certificate") or {}).get("sha256"),
            "reverse_of": record.get("reverse_of"),
            "namespaces": NAMESPACES, "slot_order_rule": SLOT_ORDER_RULE,
            "tolerance_log10": LOCALITY_TOLERANCE_LOG10}
    if check_upstream:
        why = upstream_refusal(record, record_certificate, by_id)
        if why:
            return finish({**head, "status": "refused",
                           "refusal": {"by": "upstream", "status": status_of(why),
                                       "reason": why},
                           "templates": [], "contexts": []})
    if not members:
        return finish({**head, "status": "refused",
                       "refusal": {"by": "upstream", "status": "refused:unmapped",
                                   "reason": "no provider member carries the record's channel"},
                       "templates": [], "contexts": []})
    for case, row, _ctx in members:
        if row.get("status") != "fixed_point":
            return finish({**head, "status": "refused",
                           "refusal": {"by": "upstream", "status": f"refused:window",
                                       "reason": f"member {case.label[-80:]}: {row.get('status')}"
                                                 f" {str(row.get('reason', ''))[:160]}",
                                       "window_status": row.get("status"),
                                       "covered": row.get("covered")},
                           "templates": [], "contexts": []})
    rec_r = record_graph_rows(record, "reactant_graphs")
    rec_species = counter.species_of_rows(record_graph_rows(record, "product_graphs"))
    reverse = None
    cands, ctx_rows = [], []
    for case, row, (results, merges) in members:
        case_species = counter.species_of_rows(dict(case.witness.reactant_rows()))
        rec_reac = counter.species_of_rows(rec_r)
        rev = case_species != rec_reac
        if rev and counter.species_of_rows(dict(case.witness.product_rows())) != rec_reac:
            return finish({**head, "status": "refused",
                           "refusal": {"by": "emission", "status": "refused:orientation",
                                       "reason": "the member's witness carries neither the "
                                                 "record's reactant nor its product species"},
                           "templates": [], "contexts": []})
        reverse = rev
        c, rows = context_candidates(case, results, merges, family=fam.label, pathway=pathway,
                                     direction=direction, reverse=rev, providers=providers)
        cands.extend(c)
        ctx_rows.extend(dict(r, member=case.label[-120:]) for r in rows)
    cands = dedupe(cands)
    cpath = count_pathway(fam, pathway, reverse)
    # item 2 on the record's own witness: the base templates jointly
    base = [c for c in cands if c.base and c.refusal is None and c.cut is not None]
    if not base:
        return finish({**head, "status": "refused",
                       "refusal": {"by": "emission", "status": "refused:no_base_template",
                                   "reason": "; ".join(sorted({c.refusal or "" for c in cands
                                                               if c.base}))[:300]},
                       "templates": [], "contexts": ctx_rows})
    # the record's own witness is a proxy: a base template whose window or
    # port exteriors reach into the padding beyond the proxy's artificial
    # ends cannot apply to it (the record row is then not applicable and
    # each template's weight comes from its emitted witness)
    beyond = [c for c in base if reaches_beyond(c.cut, c.witness, c.unextended)]
    record_count, weights, cal, record_map = None, None, None, {}
    if not beyond:
        probes = [c.cut.matcher(1.0, f"probe{k}") for k, c in enumerate(base)]
        count, inside, _ch = count_on(probes, rec_r, None, family=fam, direction=direction,
                                      pathway=cpath, counter=counter,
                                      product_species=rec_species, target=False)
        record_count = count.as_dict()
        if count.matcher == 0:
            raise EmissionError(f"T_ch = 0 on the record's own witness: the base template(s) "
                                f"make no transition of the record's channel (rmg_moves: "
                                f"{count.rmg})")
        if not count.equal:
            return finish({**head, "status": "refused",
                           "refusal": {"by": "emission", "status": "refused:count_mismatch",
                                       "reason": f"record witness: matcher {count.matcher} vs "
                                                 f"rmg_moves {count.rmg} transitions",
                                       "matcher": count.matcher, "rmg": count.rmg},
                           "counts": {"record": record_count}, "templates": [],
                           "contexts": ctx_rows})
        weights = weight_table(rate["T"], rate["k"], rate["ssa_multiplier"], count.matcher)
        cal = calibration(rate["T"], rate["k"], rate["ssa_multiplier"], weights,
                          count.matcher)
        for tr, _key in inside:
            record_map.setdefault(tr.template, [list(tr.mapping)])
    templates = []
    for k, cand in enumerate(cands):
        templates.append(emit_candidate(
            cand, record=record, family=fam, direction=direction, pathway=pathway,
            cpath=cpath, reverse=reverse, counter=counter, rates=rates, weights=weights,
            temperatures=rate["T"], providers=providers, ceiling=ceiling,
            owner_rate_source=owner_rate_source, record_map=next((record_map.get(f"probe{i}") for i, c in enumerate(base)
                             if c is cand), None)))
    status, refusal = record_status(templates)
    out = {**head, "status": status, "refusal": refusal, "member_reverse": bool(reverse),
           "counts": {"record": record_count},
           "weight_source": "record" if weights is not None else "emitted_witness",
           "weight": {"T": rate["T"], "w": weights} if weights is not None else None,
           "calibration": cal, "templates": templates, "contexts": ctx_rows}
    return finish(out)


def record_status(templates) -> tuple:
    """(status, refusal) of a record from its templates: ``emitted`` when the
    base template carries the rate (emitted, or merged into an equal template
    of another record) and every template does; ``partial`` when only the
    base does; else ``refused`` naming the base templates' refusals."""
    live = ("emitted", "merged")
    base_ok = any(t["status"] in live and t["base"] for t in templates)
    if not base_ok:
        reasons = [t["refusal"]["status"] for t in templates if t["base"] and t.get("refusal")]
        named = next((x for x in ("refused:collision_grid", "refused:collision", "refused:cycle",
                                   "refused:non_transferable") if x in reasons), None)
        collision = named is not None and named.startswith("refused:collision")
        return "refused", {"by": "collision" if collision else "emission",
                           "status": named or "refused:base_template",
                           "reason": "; ".join(reasons)[:300]}
    if all(t["status"] in live for t in templates):
        return "emitted", None
    return "partial", None


def emit_candidate(cand, *, record, family, direction, pathway, cpath, reverse, counter,
                   rates, weights, temperatures, providers=None, ceiling=None,
                   record_map=None, owner_rate_source=None) -> dict:
    """Count (item 2, on the context's own witness), locality (item 3) and
    the payload of one candidate."""
    entry = {"member": cand.member[-120:], "context_index": cand.context_index,
             "context": [[p[0], p[1], e] for p, e in cand.assignment if e],
             "base": cand.base, "sources": [[m[-120:], i] for m, i in cand.sources],
             "certificate": cand.certificate}
    if cand.cut is None or cand.refusal:
        entry.update(status="refused", refusal={"status": status_of(cand.refusal),
                                                "reason": cand.refusal})
        return entry
    cut = cand.cut
    entry.update(template_payload(cut, cand.witness, cand.window_atoms))
    if record_map:
        entry["maps"]["record"] = record_map
    template = cut.matcher(1.0, "t")
    rows = []
    ok = True
    reasons = []
    wset = locality_witnesses(cand.witness, providers, ceiling=ceiling)
    beyond = reaches_beyond(cut, cand.witness, cand.unextended)
    jobs = [("context", (), cand.unextended)] + [(x.label, x.vector, x) for x in wset]
    for label, vector, item in jobs:
        if label == "context" and beyond:
            rows.append({"label": label, "vector": [], "status": "not_applicable",
                         "reason": "the window or a port exterior reaches the padding beyond "
                                   "the unextended witness's artificial ends"})
            continue
        if isinstance(item, LocalityWitness):
            if item.status != "ok":
                rows.append({"label": label, "vector": list(vector), "status": item.status,
                             "reason": item.reason[:200]})
                ok = False
                reasons.append(f"{label}: {item.status}")
                continue
            witness = item.witness
        else:
            witness = item
        row = witness_row(template, witness, label=label, vector=vector, family=family,
                          direction=direction, pathway=cpath, reverse=reverse, counter=counter,
                          rates=rates, record=record, owner_rate_source=owner_rate_source)
        rows.append(row)
        if row["status"] != "ok":
            ok = False
            reasons.append(f"{label}: {row['status']}: {row.get('reason', '')[:120]}")
    emitted_row = next((r for r in rows if r["label"] == "emitted"), None)
    if weights is not None:
        record_row = {"label": "record", "vector": [], "status": "ok", "w": list(weights),
                      "source": "record k_table / T_R"}
        weight_source, cal = "record", None
    else:
        record_row = {"label": "record", "vector": [], "status": "not_applicable",
                      "reason": "a base template reaches beyond the record's proxy graph"}
        weight_source, cal = "emitted_witness", None
        if emitted_row is not None and emitted_row["status"] == "ok":
            weights = list(emitted_row["w"])
            temperatures = list(emitted_row["T"])
            cal = calibration(temperatures, emitted_row["k"], emitted_row["ssa"], weights,
                              emitted_row["counts"]["matcher"])
    good = [r for r in rows if r["status"] == "ok"] + \
        ([record_row] if record_row["status"] == "ok" else [])
    per, worst = spread(good, temperatures)
    entry["counts"] = {"context": next((r.get("counts") for r in rows
                                        if r["label"] == "context"), None),
                       "emitted": next((r.get("counts") for r in rows
                                        if r["label"] == "emitted"), None)}
    entry["weight_source"] = weight_source
    entry["weight"] = {"T": list(temperatures), "w": list(weights)} if weights else None
    if cal is not None:
        entry["calibration"] = cal
    entry["locality"] = {"rows": [record_row] + rows, "spread_log10": per,
                         "max_spread_log10": worst, "tolerance_log10": LOCALITY_TOLERANCE_LOG10,
                         "accepted": ok and worst <= LOCALITY_TOLERANCE_LOG10}
    if not ok:
        mismatch = [r for r in rows if r["status"] == "refused:count_mismatch"]
        if mismatch:
            r0 = mismatch[0]
            entry.update(status="refused", refusal={
                "status": "refused:count_mismatch",
                "reason": f"{r0['label']}: matcher {r0['counts']['matcher']} vs rmg_moves "
                          f"{r0['counts']['rmg']}", "matcher": r0["counts"]["matcher"],
                "rmg": r0["counts"]["rmg"]})
        else:
            entry.update(status="refused", refusal={"status": "refused:locality_witness",
                                                    "reason": "; ".join(reasons)[:400]})
    elif any((r.get("orbit") or {}).get("problem") for r in rows):
        bad = next(r for r in rows if (r.get("orbit") or {}).get("problem"))
        entry.update(status="refused", refusal={
            "status": "refused:non_transferable",
            "reason": f"{bad['label']}: {bad['orbit']['problem']}"[:400]})
    elif worst > LOCALITY_TOLERANCE_LOG10:
        entry.update(status="refused", refusal={
            "status": "refused:non_local_kinetics",
            "reason": f"log10 w spread {worst:.4g} > {LOCALITY_TOLERANCE_LOG10}"})
    else:
        entry.update(status="emitted", refusal=None)
    return entry


def witness_row(template, witness, *, label, vector, family, direction, pathway, reverse,
                counter, rates, record, owner_rate_source=None) -> dict:
    """One locality row: counts by both counters on the witness's own target,
    production k, ``w``."""
    r, p = template_sides(witness, reverse)
    try:
        count, inside, _ch = count_on([template], r, p, family=family, direction=direction,
                                      pathway=pathway, counter=counter)
    except wm.WindowMatcherError as error:
        return {"label": label, "vector": list(vector), "status": "refused:matcher",
                "reason": str(error)[:200]}
    row = {"label": label, "vector": list(vector), "counts": count.as_dict(),
           "heavy_atoms": witness.heavy_count()}
    if count.matcher == 0:
        raise EmissionError(f"T_ch = 0 on witness {label!r}: the template makes no transition "
                            f"of the witness's own channel (rmg_moves: {count.rmg})")
    if not count.equal:
        row.update(status="refused:count_mismatch",
                   reason=f"matcher {count.matcher} vs rmg_moves {count.rmg}")
        return row
    if not count.target_found:
        row.update(status="refused:target", reason="the witness's own mapped product is not "
                                                   "among the template's transitions")
        return row
    rs = (record or {}).get("rate_source")
    result = rates(RateRequest(family.label if hasattr(family, "label") else family,
                               (record or {}).get("template", ""),
                               (record or {}).get("orientation", "as_generated"),
                               dict(witness.reactant_rows()), dict(witness.product_rows()),
                               witness, label, tuple(vector), rs, owner_rate_source))
    if result.problems:
        row.update(status="refused:rebind", reason="; ".join(result.problems)[:300])
        return row
    w = weight_table(result.temperatures, result.k, result.ssa_multiplier, count.matcher)
    if True:
        # transferability: the channel's T_ch transitions share one weight, so
        # they must be one orbit of the reactants' automorphism group
        row["orbit"] = {"transitions": len(inside),
                        "problem": orbit_problem(r, [tr.edits for tr, _k in inside])}
    row.update(status="ok", T=list(result.temperatures), k=list(result.k),
               ssa=result.ssa_multiplier, w=w,
               provenance=result.provenance)
    return row


# ---------------------------------------------------------------------------
# items 6, 7: payload, validation, projection, record attachment
# ---------------------------------------------------------------------------
def template_payload(cut: CutTemplate, witness=None, window_atoms=None) -> dict:
    """The serialized template: executable content (hashed), metadata, maps."""
    atoms = []
    for i, a in enumerate(cut.atoms):
        boundary, labels, wport = cut.meta[i]
        atoms.append(dict(_atom_json(a), index=i, boundary=boundary, labels=list(labels),
                          witness_port=list(wport) if wport else None))
    edits = {"atoms": [[i, {k: v for k, v in a[8]}] for i, a in enumerate(cut.atoms) if a[8]],
             "bonds": [[i, j, after] for i, j, before, after in cut.bonds if before != after]}
    participants = [{"index": p, "atoms": [i for i, a in enumerate(cut.atoms) if a[0] == p],
                     "roles": sorted({lab for i, a in enumerate(cut.atoms) if a[0] == p
                                      for lab in cut.meta[i][1]})}
                    for p in range(cut.participants)]
    ports = [{"atom": i, "exterior": [list(x) for x in cut.atoms[i][7]],
              "label": cut.meta[i][0]} for i in cut.ports]
    window, maps = None, {"witness": list(cut.witness_atoms), "components": list(cut.components)}
    if witness is not None:
        window = {"atoms": len(window_atoms or ()),
                  "depths": [[list(k), d] for k, d in witness.extension_counts]}
        maps["window_id"] = hashlib.sha256(canonical_json(
            [witness.label, sorted(window_atoms or ())]).encode()).hexdigest()
        maps["witness_label"] = witness.label[:300]
    return {"name": cut.name, "key": cut.key, "family": cut.family,
            "direction": cut.direction, "pathway": list(cut.pathway),
            "executable": cut.executable(),
            "atoms": atoms, "bonds": [[i, j, x, y] for i, j, x, y in cut.bonds],
            "edits": edits, "participants": participants, "ports": ports, "predicates": [],
            "whole_species": cut.whole_species, "maps": maps, "window": window}


def finish(payload: dict) -> dict:
    """Seal a payload: sha256 over everything but the maps and the hash."""
    payload = json.loads(canonical_json(payload))
    payload["sha256"] = payload_sha256(payload)
    return payload


def _strip(payload):
    out = {k: v for k, v in payload.items() if k != "sha256"}
    out["templates"] = [{k: v for k, v in t.items() if k != "maps"}
                        for t in out.get("templates", [])]
    return out


def payload_sha256(payload: dict) -> str:
    return sha256_json(_strip(payload))


def validate_payload(payload: dict) -> None:
    """Shape and integrity of an ``executable_template`` payload (raises
    ``EmissionError``)."""
    if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
        raise EmissionError("executable_template has an unknown schema")
    if payload.get("status") not in RECORD_STATUSES:
        raise EmissionError(f"executable_template status {payload.get('status')!r}")
    if payload.get("sha256") != payload_sha256(payload):
        raise EmissionError("executable_template sha256 does not match its content")
    if payload["status"] == "refused" and not (payload.get("refusal") or {}).get("status"):
        raise EmissionError("a refused executable_template needs a refusal status")
    for t in payload.get("templates", []):
        if t.get("status") not in TEMPLATE_STATUSES:
            raise EmissionError(f"template status {t.get('status')!r}")
        if t["status"] == "refused":
            if not (t.get("refusal") or {}).get("status"):
                raise EmissionError("a refused template needs a refusal status")
            continue
        if t["status"] == "merged":
            into = t.get("merged_into") or {}
            if not into.get("record") or into.get("key") != t.get("key") or \
                    sha256_json(t["executable"]) != t["key"]:
                raise EmissionError("a merged template needs its key and owning record")
            continue
        ex = t["executable"]
        if sha256_json(ex) != t["key"]:
            raise EmissionError("template key does not match its executable content")
        n = len(ex["atoms"])
        for i, j, _x, _y in ex["bonds"]:
            if not (0 <= i < j < n):
                raise EmissionError("template bond indexes out of range")
        edited = {i for i, a in enumerate(ex["atoms"]) if a["final"]} | {
            x for i, j, before, after in ex["bonds"] if before != after for x in (i, j)}
        if any(ex["atoms"][i]["kind"] == "port" for i in edited):
            raise EmissionError("an edit touches a port")
        w = t["weight"]["w"]
        if not all(math.isfinite(x) and x > 0 for x in w):
            raise EmissionError("non-finite or non-positive weight")
        loc = t["locality"]
        if not loc["accepted"] or loc["max_spread_log10"] > LOCALITY_TOLERANCE_LOG10:
            raise EmissionError("an emitted template without an accepted locality table")
        counts = (t.get("counts") or {}).get("context") or {}
        if counts and counts.get("matcher") != counts.get("rmg"):
            raise EmissionError("an emitted template with unequal counts")
        matcher_from_payload(t, t["weight"]["T"][0])


def matcher_from_payload(template: dict, temperature: float, *, name=None) -> wm.Template:
    """The temperature-specific projection of one emitted template to a
    ``window_matcher.Template`` (``k`` interpolated linearly in ln w on the
    grid; extrapolation refused)."""
    ex = template["executable"]
    ts, ws = template["weight"]["T"], template["weight"]["w"]
    weight = _interpolate(ts, ws, float(temperature))
    atoms = [wm.WindowAtom(element=a["element"], participant=a["participant"], kind=a["kind"],
                           exterior=tuple(tuple(x) for x in a["exterior"]),
                           isotope=a["isotope"], charge=a["charge"], radical=a["radical"],
                           lone_pairs=a["lone_pairs"]) for a in ex["atoms"]]
    bonds = [(i, j, x) for i, j, x, _y in ex["bonds"] if x]
    edits = [wm.BondEdit(i, j, y) for i, j, x, y in ex["bonds"] if x != y]
    edits += [wm.AtomEdit(i, **a["final"]) for i, a in enumerate(ex["atoms"]) if a["final"]]
    return wm.Template(name or template["name"], template.get("family", "") or "", atoms,
                       bonds, edits, weight=weight)


def _interpolate(ts, ws, t):
    ts = [float(x) for x in ts]
    if not ts or t < ts[0] or t > ts[-1]:
        raise EmissionError(f"temperature {t} K outside the weight table (extrapolation refused)")
    for i in range(len(ts) - 1):
        if ts[i] <= t <= ts[i + 1]:
            if ts[i + 1] == ts[i]:
                return ws[i]
            f = (t - ts[i]) / (ts[i + 1] - ts[i])
            return math.exp(math.log(ws[i]) * (1 - f) + math.log(ws[i + 1]) * f)
    return ws[-1]


def ownership_problem(record, payload: dict) -> Optional[str]:
    """Why ``payload`` is not this record's (None: it is): a payload is
    emitted for one record, named by its event id, family and orientation."""
    get = record.get if isinstance(record, dict) else (lambda k: getattr(record, k, None))
    for field_, mine in (("record_event_id", get("event_id")), ("family", get("family")),
                         ("orientation", get("orientation"))):
        if payload.get(field_) != mine:
            return (f"the executable_template belongs to {field_} {payload.get(field_)!r}, "
                    f"not this record's {mine!r}")
    return None


def attach(record, payload: dict):
    """An ``EventRecord`` with ``executable_template`` set (schema
    ``rewrite_record/0.2``; the event id is unchanged: the template is not
    part of it).  Only the record's own payload attaches
    (``ownership_problem``)."""
    from dataclasses import replace
    from rmgpy.kmc.event_record import EventRecord
    validate_payload(payload)
    rec = record if isinstance(record, EventRecord) else EventRecord.from_dict(dict(record))
    why = ownership_problem(rec, payload)
    if why:
        raise EmissionError(why)
    return replace(rec, schema_version=RECORD_SCHEMA, executable_template=payload)


# ---------------------------------------------------------------------------
# catalog-level consistency (I-147)
# ---------------------------------------------------------------------------
class ClassOverlapError(EmissionError):
    """One family emitted two different templates that match one site with
    one edit: the transition would be counted twice."""


class ReciprocalSplitError(EmissionError):
    """A merge that would split a reciprocal (``reverse_of``) pair of records."""


CYCLE_RTOL = 1e-12
COLLISION_STATUS = "refused:collision"


def orbit_problem(rows, edit_sets) -> Optional[str]:
    """Why the transitions (``edit_sets``: each one's ``window_matcher``
    edits) are not one orbit under the automorphism group of the witness
    reactants ``rows`` (None: they are one orbit).

    Two transitions are equivalent when an automorphism of the labelled
    reactant graph (element, isotope, charge, radicals, lone pairs, bond
    orders; the two copies of an A + A pair and their exchange included)
    maps the edits of one onto the edits of the other: the graph with each
    edit written on its atoms and bonds (what the atom ends as, the bond
    order it ends at) is isomorphic to the other's.  An orbit is transitive,
    so every transition is compared with the first."""
    edit_sets = list(edit_sets)
    if len(edit_sets) < 2:
        return None
    import networkx as nx
    from networkx.algorithms import isomorphism as iso
    from rmgpy.kmc import melt_state as ms
    base = nx.Graph()
    for a, row in rows.items():
        base.add_node(a, lab=repr(tuple(row[:5])))
    for a, row in rows.items():
        for n, o in row[5]:
            if a < n:
                base.add_edge(a, n, order=float(o), elab=repr((float(o), None)))

    def annotated(edits):
        g = base.copy()
        marks = defaultdict(list)
        for e in edits:
            if isinstance(e, ms.SetBond):
                before = g.edges[e.a, e.b]["order"] if g.has_edge(e.a, e.b) else 0.0
                g.add_edge(e.a, e.b, order=before, elab=repr((before, float(e.new))))
            else:
                marks[e.atom].append((type(e).__name__, e.new))
        for a, m in marks.items():
            g.nodes[a]["lab"] = g.nodes[a]["lab"] + repr(tuple(sorted(m)))
        return g

    first = annotated(edit_sets[0])
    h0 = nx.weisfeiler_lehman_graph_hash(first, node_attr="lab", edge_attr="elab", iterations=5)
    for k, edits in enumerate(edit_sets[1:], 1):
        other = annotated(edits)
        same = nx.weisfeiler_lehman_graph_hash(other, node_attr="lab", edge_attr="elab",
                                               iterations=5) == h0 and \
            iso.GraphMatcher(first, other,
                             node_match=lambda x, y: x["lab"] == y["lab"],
                             edge_match=lambda x, y: x["elab"] == y["elab"]).is_isomorphic()
        if not same:
            return (f"the channel's transitions are not one orbit of the witness's "
                    f"automorphism group: transition 0 and transition {k} are not related "
                    f"by an automorphism of the reactants")
    return None


def _grid_problem(a, b) -> Optional[str]:
    a, b = list(a), list(b)
    if a != b:
        return f"temperature grids differ ({len(a)} points {a[:1]}..{a[-1:]} vs {len(b)} points " \
               f"{b[:1]}..{b[-1:]})"
    return None


def _weight_problem(payload) -> Optional[str]:
    """A non-finite or non-positive weight, or an empty channel, in a
    payload's emitted template (compile errors)."""
    for t in payload.get("templates", []):
        if t.get("status") != "emitted":
            continue
        w = (t.get("weight") or {}).get("w") or []
        grid = (t.get("weight") or {}).get("T") or []
        if not w or len(w) != len(grid):
            return f"template {t.get('key', '')[:12]}: weight table is empty or ragged"
        for T, x in zip(grid, w):
            if not isinstance(x, (int, float)) or not math.isfinite(x) or x <= 0:
                return (f"template {t.get('key', '')[:12]}: non-finite or non-positive weight "
                        f"{x!r} at {T} K")
        counts = (t.get("counts") or {}).get("context") or {}
        if counts and counts.get("matcher") == 0:
            return f"template {t.get('key', '')[:12]}: T_ch = 0 on its context witness"
    counts = (payload.get("counts") or {}).get("record") or {}
    if counts and counts.get("matcher") == 0:
        return "T_ch = 0 on the record's own witness"
    return None


def _live_bases(payload, index=None):
    """The record's emitted base templates; with ``index`` ((record, key) ->
    template) a merged base template stands as its owner's structure at the
    copy's own weight (the cycle checks the records as compiled)."""
    out = []
    for t in payload.get("templates", []):
        if not t.get("base"):
            continue
        if t.get("status") == "emitted":
            out.append(t)
        elif t.get("status") == "merged" and index is not None:
            owner = index.get((t["merged_into"]["record"], t["merged_into"]["key"]))
            if owner is not None and owner.get("status") == "emitted":
                # the owner's structure with the copy's own (calibrated) weight
                out.append(dict(owner, weight=t["weight"]))
    return out


def _reciprocal_partitions(out, groups):
    """Check that the partition of templates by key is the same on both sides
    of every reciprocal record pair; return {(record, key): (reverse record,
    reverse key)} for the base templates with exactly one live partner.

    ``groups``: key -> [(record, payload index, template index)] of every
    emitted template.  Only base templates take part, and only records whose
    reverse record has exactly one live base template."""
    by_record = {p.get("record_event_id"): p for p in out}
    links, base_members = {}, defaultdict(set)
    for key, members in groups.items():
        for rid, i, j in members:
            if out[i]["templates"][j].get("base"):
                base_members[key].add(rid)
    for rid, p in by_record.items():
        mine = _live_bases(p)
        rev = p.get("reverse_of")
        partner = by_record.get(rev)
        if len(mine) != 1 or partner is None:
            continue
        theirs = _live_bases(partner)
        if len(theirs) == 1:
            links[(rid, mine[0]["key"])] = (rev, theirs[0]["key"])
    for key in sorted(base_members):
        mine = {r for r in base_members[key] if (r, key) in links}
        targets = {links[(r, key)][1] for r in mine}
        if len(targets) > 1:
            raise ReciprocalSplitError(
                f"reciprocal pairs split: the records {sorted(mine)} share template "
                f"{key[:16]} but their reverse records fall into {len(targets)} templates "
                f"{sorted(t[:16] for t in targets)}")
        for target in targets:
            partners = {links[(r, key)][0] for r in mine}
            others = {r for r in base_members[target] if (r, target) in links}
            if partners != others:
                raise ReciprocalSplitError(
                    f"reciprocal pairs split: template {key[:16]} (records {sorted(mine)}) "
                    f"pairs with template {target[:16]} (records {sorted(others)}), whose "
                    f"records are not exactly the reverses of the first set")
    return links


def resolve_collisions(payloads, tolerance_log10: float = LOCALITY_TOLERANCE_LOG10) -> tuple:
    """One catalog's payloads with equal executable keys made one template.

    Across every emitted template of every payload (one per record), a key
    emitted more than once would fire every matched transition once per
    copy.  Copies on one temperature grid whose log10 weights have a
    DIAMETER (max - min over all the copies, at each grid temperature) of at
    most ``tolerance_log10`` are merged, whatever the input order: the copy of
    the least record event id stays ``emitted`` at its own weight (the
    representative; the weights of the others differ by less than the
    locality tolerance) and lists every source record (``collision``); the
    others become ``merged`` (``merged_into``).  Copies beyond the diameter
    are all refused ``refused:collision`` (keys, records, weights); copies on
    different grids ``refused:collision_grid``.

    Reciprocal pairs (``reverse_of``) merge together or not at all
    (``ReciprocalSplitError``); the emitted base templates of a linked pair
    carry ``reverse`` (the reverse template key and records), bijectively.
    A non-finite or non-positive weight, or ``T_ch = 0``, in any emitted
    template is a compile error (``EmissionError``).  Record statuses follow
    (``record_status``).  Returns (new payloads, collision table: one row per
    key with more than one copy)."""
    for p in payloads:
        why = _weight_problem(p)
        if why:
            raise EmissionError(f"record {str(p.get('record_event_id'))[:16]}: {why}")
    out = [json.loads(canonical_json(p)) for p in payloads]
    groups = defaultdict(list)
    for i, p in enumerate(out):
        for j, t in enumerate(p.get("templates", [])):
            if t.get("status") == "emitted":
                groups[t["key"]].append((p.get("record_event_id") or "", i, j))
    links = _reciprocal_partitions(out, groups)
    table = []
    touched = set()
    for key in sorted(groups):
        copies = sorted(groups[key])
        if len(copies) < 2:
            continue
        weights = [out[i]["templates"][j]["weight"] for _r, i, j in copies]
        records = [r for r, _i, _j in copies]
        grid = next((_grid_problem(weights[0]["T"], w["T"]) for w in weights[1:]
                     if _grid_problem(weights[0]["T"], w["T"])), None)
        diameter = None
        if grid is None:
            logs = [[math.log10(x) for x in w["w"]] for w in weights]
            diameter = max(max(c) - min(c) for c in zip(*logs))
            spread = max((max(c) - min(c)) / max(c) for c in zip(*[w["w"] for w in weights]))
            merge = diameter <= tolerance_log10
            disposition = "merged" if merge else COLLISION_STATUS
        else:
            spread, merge, disposition = None, False, COLLISION_STATUS + "_grid"
        row = {"key": key, "records": records, "copies": len(copies),
               "log10_diameter": diameter, "rel_spread": spread, "disposition": disposition,
               "tolerance_log10": tolerance_log10,
               "w_700": [w["w"][list(w["T"]).index(700.0)] if 700.0 in list(w["T"]) else None
                         for w in weights]}
        table.append(row)
        owner_record, oi, oj = copies[0]
        for r, i, j in copies:
            t = out[i]["templates"][j]
            touched.add(i)
            if merge and (i, j) == (oi, oj):
                t["collision"] = {"merged": True, "records": records,
                                  "log10_diameter": diameter, "rel_spread": spread,
                                  "tolerance_log10": tolerance_log10,
                                  "weight_rule": "the least record event id's weight"}
            elif merge:
                t["status"] = "merged"
                t["merged_into"] = {"record": owner_record, "key": key}
            else:
                t["status"] = "refused"
                reason = (f"{len(copies)} copies of key {key[:16]} on different grids: {grid}"
                          if grid else
                          f"{len(copies)} copies of key {key[:16]} with a log10 weight diameter "
                          f"{diameter:.4g} > {tolerance_log10}")
                t["refusal"] = {"status": disposition, "reason": reason, "keys": [key],
                                "records": records,
                                "weights": [{"record": rr, "T": w["T"], "w": w["w"]}
                                            for rr, w in zip(records, weights)]}
    # the bijective reverse link of every linked base template
    live = defaultdict(set)
    for i, p in enumerate(out):
        for t in p.get("templates", []):
            if t.get("status") in ("emitted", "merged") and t.get("base"):
                rid = p.get("record_event_id")
                if t["status"] == "merged":
                    live[t["merged_into"]["key"]].add(rid)
                else:
                    live[t["key"]].add(rid)
    for i, p in enumerate(out):
        rid = p.get("record_event_id")
        for t in p.get("templates", []):
            link = links.get((rid, t.get("key")))
            if link is None or t.get("status") not in ("emitted", "merged") or not t.get("base"):
                continue
            partners = sorted({links[(r, t["key"])][0] for r in live[t["key"]]
                               if (r, t["key"]) in links}) if t["status"] == "emitted" \
                else [link[0]]
            t["reverse"] = {"key": link[1], "records": partners}
            touched.add(i)
    for i in touched:
        out[i]["status"], out[i]["refusal"] = record_status(out[i]["templates"])
        out[i] = finish(out[i])
    return out, table


def check_cycles(payloads, records, *, counter, rtol: float = CYCLE_RTOL) -> tuple:
    """The equilibrium cycle of every reciprocal pair on its witness.

    For records F and R linked by ``reverse_of`` (both payloads with live
    base templates that apply to the records' own witnesses), on F's reactant
    witness ``S_F`` the transitions of F's base templates in F's channel
    (``T_F`` of them) and on R's reactant witness ``S_R`` those of R's
    (``T_R``):

    * bijection: the inverse edit of each F transition, applied on its
      product state, is a transition of R's templates (exactly one returns
      ``S_F``), and the inverse edit of each R transition on its product state
      is a transition of F's templates (exactly one returns ``S_R``);
    * ratio: ``(sum of the F transitions' template weights) / (sum of the R
      transitions' weights)`` equals ``(ssa_F k_F) / (ssa_R k_R)`` from the
      two records' compiled ``k_table`` (the class Kc as compiled, with the
      records' own multipliers) to ``rtol`` at every grid temperature.

    A failure refuses every emitted template of both records
    (``refused:cycle``, with the numbers).  Pairs whose templates reach the
    padding beyond a record's proxy graph (no ``counts.record``) are
    ``not_applicable``; a pair whose grids differ is refused by name.
    Returns (new payloads, one table row per pair)."""
    out = [json.loads(canonical_json(p)) for p in payloads]
    by_record = {p.get("record_event_id"): i for i, p in enumerate(out)}
    index = {(p.get("record_event_id"), t.get("key")): t for p in out
             for t in p.get("templates", []) if t.get("status") == "emitted"}
    table, refuse = [], {}
    for rid in sorted(by_record):
        p = out[by_record[rid]]
        rev = p.get("reverse_of")
        if rev is None or rev not in by_record or rid > rev:
            continue
        q = out[by_record[rev]]
        row = {"forward": rid, "reverse": rev}
        table.append(row)
        if q.get("reverse_of") != rid:
            row.update(status="refused:cycle", reason="the reverse link is not reciprocal")
            refuse[rid] = refuse[rev] = row["reason"]
            continue
        if not _live_bases(p, index) or not _live_bases(q, index):
            row.update(status="skipped", reason="a record has no live base template")
            continue
        if (p.get("counts") or {}).get("record") is None or \
                (q.get("counts") or {}).get("record") is None or \
                rid not in records or rev not in records:
            row.update(status="not_applicable",
                       reason="a base template reaches beyond a record's proxy graph")
            continue
        try:
            result = _cycle(p, q, records[rid], records[rev], counter, rtol, index)
        except EmissionError as error:
            result = {"status": "refused:cycle", "reason": str(error)}
        row.update(result)
        if row["status"] != "ok":
            refuse[rid] = refuse[rev] = f"pair {rid[:16]}/{rev[:16]}: {row['reason']}"
    dead = {}
    for rid, reason in refuse.items():
        for t in out[by_record[rid]].get("templates", []):
            if t.get("status") == "emitted":
                dead[(rid, t["key"])] = reason
    changed = set()
    for i, p in enumerate(out):
        rid = p.get("record_event_id")
        for t in p.get("templates", []):
            if t.get("status") == "emitted" and (rid, t["key"]) in dead:
                why = dead[(rid, t["key"])]
            elif t.get("status") == "merged" and (
                    t["merged_into"]["record"], t["merged_into"]["key"]) in dead:
                why = dead[(t["merged_into"]["record"], t["merged_into"]["key"])]
            else:
                continue
            t["status"] = "refused"
            t["refusal"] = {"status": "refused:cycle", "reason": why[:600]}
            t.pop("merged_into", None)
            changed.add(i)
    for i in sorted(changed):
        out[i]["status"], out[i]["refusal"] = record_status(out[i]["templates"])
        out[i] = finish(out[i])
    return out, table


def _record_channel(payload, record, counter, index):
    """(rows, state, [(transition, template payload)]): the record's reactant
    witness and the transitions of the payload's live base templates in the
    record's channel."""
    fam = counter.family(payload["family"])
    pathway = tuple(payload.get("pathway") or ())
    rows = record_graph_rows(record, "reactant_graphs")
    species = counter.species_of_rows(record_graph_rows(record, "product_graphs"))
    bases = _live_bases(payload, index)
    templates = [matcher_from_payload(t, t["weight"]["T"][0], name=f"t{k}")
                 for k, t in enumerate(bases)]
    channel = counter.channel(fam, molecules_of(rows), payload["direction"],
                              count_pathway(fam, pathway, bool(payload.get("member_reverse"))),
                              species)
    state = melt_state(rows)
    keyset = set(channel.keys)
    inside = [(tr, bases[int(tr.template[1:])])
              for tr in wm.match(templates, state).transitions
              if rm.product_key(wm.product_rows(state, tr)) in keyset]
    return rows, state, inside


def _norm_row(row):
    return (row[0], row[1], int(row[2]), int(row[3]), int(row[4]),
            tuple(sorted((int(n), float(o)) for n, o in row[5])))


def _returns(rows, state, inside, payload, index):
    """Per transition of ``inside`` (on ``state``/``rows``): the number of
    transitions of ``payload``'s live base templates on its product state
    that give the original components back."""
    templates = [matcher_from_payload(t, t["weight"]["T"][0], name=f"t{k}")
                 for k, t in enumerate(_live_bases(payload, index))]
    first = first_form_rows(rows)
    counts = []
    for tr, _t in inside:
        touched = wm.product_rows(state, tr)
        after = dict(rows)
        after.update(touched)
        check = melt_state(after)
        back = 0
        for rt in wm.match(templates, check).transitions:
            got = wm.product_rows(check, rt)
            if set(got) == set(touched) and \
                    all(_norm_row(got[a]) == _norm_row(first[a]) for a in got):
                back += 1
        counts.append(back)
    return counts


def _cycle(p, q, rec_p, rec_q, counter, rtol, index) -> dict:
    rows_p, state_p, in_p = _record_channel(p, rec_p, counter, index)
    rows_q, state_q, in_q = _record_channel(q, rec_q, counter, index)
    if not in_p or not in_q:
        raise EmissionError("a record's base templates make no transition of its channel on "
                            "its own witness")
    back_p = _returns(rows_p, state_p, in_p, q, index)
    back_q = _returns(rows_q, state_q, in_q, p, index)
    bij = {"forward_to_reverse": all(n == 1 for n in back_p),
           "reverse_to_forward": all(n == 1 for n in back_q),
           "forward_returns": back_p, "reverse_returns": back_q}
    base = {"T_forward": len(in_p), "T_reverse": len(in_q), "bijection": bij}
    if not (bij["forward_to_reverse"] and bij["reverse_to_forward"]):
        return dict(base, status="refused:cycle", reason=(
            f"inverse edit: each forward transition must be undone by exactly one reverse "
            f"transition and vice versa (forward {back_p}, reverse {back_q})"))
    grids = {tuple(t["weight"]["T"]) for _tr, t in in_p + in_q}
    grids.add(tuple(rec_p["k_table"]["T"]))
    grids.add(tuple(rec_q["k_table"]["T"]))
    if len(grids) != 1:
        return dict(base, status="refused:cycle", reason="ratio: the pair's templates and "
                    "records are on different temperature grids")
    grid = list(next(iter(grids)))
    lhs = [math.fsum(t["weight"]["w"][i] for _tr, t in in_p) for i in range(len(grid))]
    rhs = [math.fsum(t["weight"]["w"][i] for _tr, t in in_q) for i in range(len(grid))]
    kc = [float(rec_p.get("ssa_multiplier", 1.0)) * kf / (float(rec_q.get("ssa_multiplier", 1.0)) * kr)
          for kf, kr in zip(rec_p["k_table"]["k"], rec_q["k_table"]["k"])]
    errors = [abs(a / b / c - 1.0) for a, b, c in zip(lhs, rhs, kc)]
    worst = max(errors)
    i = errors.index(worst)
    base.update(max_ratio_error=worst, worst_T=grid[i], Kc_worst=kc[i],
                k_forward_worst=rec_p["k_table"]["k"][i], k_reverse_worst=rec_q["k_table"]["k"][i],
                w_forward_sum_worst=lhs[i], w_reverse_sum_worst=rhs[i])
    if worst > rtol:
        return dict(base, status="refused:cycle", reason=(
            f"ratio: (w_f T_f)/(w_r T_r) = {lhs[i] / rhs[i]:.17g} vs Kc = {kc[i]:.17g} at "
            f"{grid[i]} K, relative error {worst:.3g} > {rtol} (T_f {len(in_p)}, T_r {len(in_q)})"))
    return dict(base, status="ok", reason="")


def class_overlaps(templates, rows_or_state) -> list:
    """Sites where two different templates of one family make one edit on
    ``rows_or_state`` (the state, or atom rows): rows of {family, templates,
    components}.  Different families making one net edit are distinct
    pathways and are not overlaps."""
    state = rows_or_state if hasattr(rows_or_state, "components") else \
        melt_state(rows_or_state)
    seen = defaultdict(set)
    for tr in wm.match(list(templates), state).transitions:
        seen[(tr.family, tr.components, tr.delta)].add(tr.template)
    return [{"family": fam, "templates": sorted(names), "components": list(comps)}
            for (fam, comps, _delta), names in sorted(seen.items(), key=lambda x: str(x[0]))
            if len(names) > 1]


def check_class_overlap(templates, rows_or_state) -> None:
    found = class_overlaps(templates, rows_or_state)
    if found:
        raise ClassOverlapError(
            f"{len(found)} site(s) where one family's templates make one edit twice, e.g. "
            f"family {found[0]['family']}: {', '.join(found[0]['templates'])}")


def loaded_templates(payloads, temperature: float, *, name=None) -> list:
    """The ``window_matcher.Template`` set a catalog loads at ``temperature``:
    every ``emitted`` template (merged copies and refused ones are not).  A
    ``merged`` template whose owning record is not among ``payloads`` would
    leave its channel without a template: refused loudly (``carry_owners``
    adds them to a subset)."""
    have = {(p.get("record_event_id"), t.get("key")) for p in payloads
            for t in p.get("templates", []) if t.get("status") == "emitted"}
    for p in payloads:
        for t in p.get("templates", []):
            if t.get("status") == "merged":
                into = t["merged_into"]
                if (into["record"], into["key"]) not in have:
                    raise EmissionError(
                        f"orphaned merged template {t['key'][:16]} of record "
                        f"{str(p.get('record_event_id'))[:16]}: its owner "
                        f"{into['record'][:16]} is not loaded")
    out = []
    for p in payloads:
        for t in p.get("templates", []):
            if t.get("status") == "emitted":
                label = (name or (lambda p, t: f"{(p.get('record_event_id') or '')[4:16]}:"
                                               f"{t['key'][:12]}"))(p, t)
                out.append(matcher_from_payload(t, temperature, name=label))
    return out


def carry_owners(selected, catalog) -> list:
    """``selected`` payloads (a subset of ``catalog``) plus, for each merged
    template in them whose owner is not selected, a payload of that one owning
    template (not the owner's other templates)."""
    have = {(p.get("record_event_id"), t.get("key")) for p in selected
            for t in p.get("templates", []) if t.get("status") == "emitted"}
    by_record = {p.get("record_event_id"): p for p in catalog}
    extra = {}
    for p in selected:
        for t in p.get("templates", []):
            if t.get("status") != "merged":
                continue
            into = t["merged_into"]
            if (into["record"], into["key"]) in have:
                continue
            owner = by_record.get(into["record"])
            mine = [x for x in (owner or {}).get("templates", [])
                    if x.get("key") == into["key"] and x.get("status") == "emitted"]
            if not mine:
                raise EmissionError(f"owner {into['record'][:16]} of merged template "
                                    f"{t['key'][:16]} is not in the catalog")
            extra[(into["record"], into["key"])] = {
                "record_event_id": into["record"], "family": owner.get("family"),
                "status": "emitted", "carried_for": p.get("record_event_id"),
                "templates": [mine[0]]}
    return list(selected) + [extra[k] for k in sorted(extra)]


# ---------------------------------------------------------------------------
# item 4: ladder
# ---------------------------------------------------------------------------
def ladder_member(case, row, *, family, counter, rates, orientation, ssa=None,
                  certifier=None) -> dict:
    """A ladder member's whole-species template (no locality test): the rate
    through the production path on the member's own species, ``T_ch`` by both
    counters on it."""
    fam = counter.family(family)
    entry = {"member": case.label[-120:], "direction": case.direction,
             "pathway": list(case.pathway), "kind": case.kind}
    if row.get("status") != "fixed_point":
        entry.update(status="refused", refusal={"by": "upstream",
                                                "status": f"refused:window",
                                                "window_status": row.get("status"),
                                                "reason": str(row.get("reason", ""))[:200]})
        return entry
    labels = labels_of(case.moves[0])
    cut = whole_species_cut(case.witness, family=fam.label, pathway=case.pathway,
                            direction=case.direction, labels=labels)
    entry.update(template_payload(cut, case.witness, set(case.witness.reactant_rows())))
    template = cut.matcher(1.0, "t")
    r, p = template_sides(case.witness)
    count, _inside, _ch = count_on([template], r, p, family=fam, direction=case.direction,
                                   pathway=tuple(case.pathway), counter=counter)
    entry["counts"] = {"own": count.as_dict()}
    if not count.equal:
        entry.update(status="refused", refusal={
            "status": "refused:count_mismatch",
            "reason": f"matcher {count.matcher} vs rmg_moves {count.rmg}",
            "matcher": count.matcher, "rmg": count.rmg})
        return entry
    result = rates(RateRequest(fam.label, ";".join(case.pathway), orientation, r, p,
                               case.witness, "ladder", ()))
    if result.problems:
        entry.update(status="refused", refusal={"status": "refused:rate",
                                                "reason": "; ".join(result.problems)[:300]})
        return entry
    w = weight_table(result.temperatures, result.k, result.ssa_multiplier, count.matcher)
    entry.update(status="emitted", refusal=None,
                 weight={"T": list(result.temperatures), "w": w},
                 rate={"k": list(result.k), "ssa_multiplier": result.ssa_multiplier,
                       "provenance": result.provenance},
                 calibration=calibration(result.temperatures, result.k, result.ssa_multiplier,
                                         w, count.matcher))
    return entry


def ladder_control(case, row, fp, *, family, counter) -> dict:
    """The L + 1 control: the member's fixed-point window cut (``fp``: its
    ``windows.FixedPointResult``) is a LOCAL template (it has a port: not
    whole-species) and covers its own event."""
    fam = counter.family(family)
    out = {"member": case.label[-120:], "status": row.get("status")}
    if row.get("status") != "fixed_point" or fp is None:
        out["verdict"] = "refused upstream"
        return out
    labels = labels_of(case.moves[0])
    cut = cut_template(fp.witness, fp.window.atoms, family=fam.label,
                       pathway=case.pathway, direction=case.direction, labels=labels)
    r, p = template_sides(case.witness)
    count, _inside, _ch = count_on([cut.matcher(1.0, "t")], r, p, family=fam,
                                   direction=case.direction, pathway=tuple(case.pathway),
                                   counter=counter)
    out.update(local=not cut.whole_species, ports=len(cut.ports), counts=count.as_dict(),
               verdict="local template covers it" if (not cut.whole_species and
                                                       count.target_found and count.equal)
               else "not covered by a local template")
    return out


# ---------------------------------------------------------------------------
# item 5 and acceptance: matcher union vs rmg_moves on a reference
# ---------------------------------------------------------------------------
def classify_reference(templates, molecules, *, family, directions, counter, owners=None):
    """Matcher union over ``templates`` (``window_matcher.Template``, unique
    names) on the reference ``molecules`` (separate components; A + A as two
    copies) against ``rmg_moves`` on the whole molecules in ``directions``
    (unimolecular families on each molecule, bimolecular on the pair).

    ``owners(move) -> (status, reason)``: the upstream/emission status of the
    record owning a kept RMG move (its channel: template and product species);
    ``BEYOND_PROXY`` explains an RMG-only move whose owner's templates reach
    beyond the reference's artificial ends.

    Returns a dict with the counts, the per-key overlap multiplicity and the
    classified one-sided keys; ``failures`` lists what fails (an unexplained
    RMG-only kept move, a matcher key RMG files only as forbidden, a missing
    required unchanged_species move, an unclassified matcher-only key)."""
    fam = counter.family(family)
    rows = {}
    for m in molecules:
        rows.update(rm.molecule_rows(m))
    state = melt_state(rows)
    result = wm.match(list(templates), state)
    matcher = defaultdict(list)
    for tr in result.transitions:
        matcher[rm.product_key(wm.product_rows(state, tr))].append(tr.template)
    rmg = defaultdict(lambda: defaultdict(list))
    n_moves = Counter()
    for direction in directions:
        two = _molecularity(fam, direction)
        groups = [[m] for m in molecules] if two == 1 else \
            [[molecules[i], molecules[j]] for i in range(len(molecules))
             for j in range(i + 1, len(molecules))]
        for group in groups:
            for mv in counter.moves(fam, group, direction):
                n_moves[mv.status] += 1
                if mv.products is None:
                    continue
                rmg[rm.product_key(mv)][mv.status].append(mv)
    kept = {k for k, v in rmg.items() if "kept" in v}
    unchanged = {k for k, v in rmg.items() if "unchanged_species" in v and "kept" not in v}
    forbidden = {k for k, v in rmg.items() if set(v) == {"forbidden"}}
    other = {k for k, v in rmg.items()} - kept - unchanged - forbidden
    mkeys = set(matcher)
    failures, rmg_only, matcher_only = [], [], []
    for k in sorted(kept - mkeys):
        mv = rmg[k]["kept"][0]
        status, reason = owners(mv) if owners else ("unowned", "no owner function")
        explained = (status.startswith("refused") and status not in EMISSION_FAILURE_STATUSES) \
            or status == BEYOND_PROXY
        rmg_only.append({"template": list(mv.template or ()), "labels": [list(x) for x in mv.labels],
                         "owner_status": status, "reason": reason[:200], "explained": explained})
        if not explained:
            failures.append(f"RMG-only kept move ({';'.join(mv.template or ())[:60]}): {status}")
    for k in sorted(mkeys - kept):
        names = sorted(set(matcher[k]))
        if k in unchanged:
            cls = "unchanged_species (matcher-only by policy)"
        elif k in forbidden:
            cls = "forbidden (policy: excluded)"
            failures.append(f"matcher produces a forbidden move ({names[0][:60]})")
        elif k in other:
            cls = "rmg filtered: " + ",".join(sorted(rmg[k]))
            failures.append(f"matcher produces a move RMG files as {sorted(rmg[k])}")
        else:
            cls = "unclassified"
            failures.append(f"unclassified matcher-only move ({names[0][:60]})")
        matcher_only.append({"templates": names, "class": cls})
    for k in sorted(unchanged - mkeys):
        mv = rmg[k]["unchanged_species"][0]
        status, reason = owners(mv) if owners else ("unowned", "")
        if not status.startswith("refused") or status in EMISSION_FAILURE_STATUSES:
            failures.append("required unchanged_species move not produced "
                            f"({';'.join(mv.template or ())[:60]})")
    for k, names in sorted(matcher.items()):
        if len(names) > 1:
            failures.append(f"one transition is produced {len(names)} times by "
                            f"{len(set(names))} loaded template(s) ({sorted(set(names))[0][:40]})")
    overlap = Counter(len(set(v)) for v in matcher.values())
    return {"rmg_status_moves": dict(n_moves),
            "rmg_kept": len(kept), "rmg_forbidden": len(forbidden),
            "rmg_unchanged_species": len(unchanged), "rmg_other": len(other),
            "matcher_transitions": sum(len(v) for v in matcher.values()),
            "matcher_keys": len(mkeys), "matched": len(kept & mkeys),
            "matcher_only": len(mkeys - kept), "rmg_only": len(kept - mkeys),
            "rmg_only_beyond_proxy": sum(1 for x in rmg_only if x["owner_status"] == BEYOND_PROXY),
            "unchanged_produced": len(unchanged & mkeys),
            "overlap_multiplicity": {str(k): v for k, v in sorted(overlap.items())},
            "overlap_keys": sum(1 for v in matcher.values() if len(set(v)) > 1),
            "rmg_directed_per_kept_key": dict(Counter(len(rmg[k]["kept"]) for k in kept)),
            "rmg_only_moves": rmg_only, "matcher_only_moves": matcher_only,
            "failures": failures}


# the owner emitted, but its templates reach beyond the proxy's artificial
# ends (``reaches_beyond``): on the capped reference they cannot match
BEYOND_PROXY = "emitted:beyond_proxy"

# emission refusals are named, but an RMG-only move they leave is not an
# upstream explanation (the acceptance fails on it)
EMISSION_FAILURE_STATUSES = ("refused:count_mismatch", "refused:non_local_kinetics",
                             "refused:locality_witness", "refused:base_template",
                             "refused:no_base_template", "refused:orientation",
                             "refused:unmapped")


def _molecularity(family, direction) -> int:
    """Reactant count of a family direction (``windows.family_constants``:
    a family whose template is a single tree root, such as R_Recombination,
    lists one template reactant for a bimolecular direction)."""
    from rmgpy.kmc import windows as win
    return dict(win.family_constants(family).directions)[direction]
