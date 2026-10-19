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
on every matched transition.  ``resolve_collisions`` merges copies whose
log10 weights have a diameter of at most ``LOCALITY_TOLERANCE_LOG10`` at
every grid temperature (the weight of the copy of least
``partition_free_identity``, never chosen by record id; every source listed) and
refuses the rest.  Reciprocal (``reverse_of``) groups are admitted or refused
together and their representatives come from one source record pair.
``check_cycles`` then tests the LOADED catalog's detailed balance on the
stored witnesses against independent RMG thermo (``ThermoKc``).
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
    {(i, j): label} with string labels).  The refinement starts from the
    executable labels alone: the metadata only breaks ties between orders of
    equal executable serialization, so equal executables get one key whatever
    their metadata.  Interchangeable hydrogen leaves (degree one, same
    neighbour and bond, same labels) are automorphic and are not branched."""
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
        if all(len(nbrs[v]) == 1 for v in target) and len({nbrs[v][0] for v in target}) == 1 \
                and len({meta_labels[v] for v in target}) == 1:
            choices = [min(target)]   # automorphic leaves: one branch serializes them all
        else:
            choices = sorted(target)
        for v in choices:
            search(refine(rank([(c, 0 if u == v else 1) for u, c in enumerate(colors)])))

    search(refine(rank(list(exec_labels))))
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

    def inverts(self, family, move, molecules, pathway) -> bool:
        """``move`` (a reading of ``molecules`` in the template direction of an
        own-reverse family) is the exact inverse of an RMG move of the record's
        recipe: some move of the same family on ``move``'s product, carrying
        the node labels ``pathway``, restores the touched reactant graph
        (atom ids kept, modulo resonance).  RMG's own reading of the record's
        reactants may sit on another tree node, and other recipes (another ring
        size) may reach the same product species; neither is the record's
        reaction."""
        family = self.family(family)
        product = move.product_rows()
        ids = set(product)
        back = {}
        for mol in molecules:
            rows = rm.molecule_rows(mol)
            if set(rows) <= ids:
                back.update(rows)
        if set(back) != ids:
            return False
        want = rm.product_key(back)
        parts = [rm.molecule_from_rows(part) for part in rm.split_rows(dict(product))]
        for mv in self.moves(family, parts, "forward"):
            if mv.status not in ("kept", "unchanged_species") or mv.products is None \
                    or mv.template is None or tuple(mv.template) != tuple(pathway):
                continue
            if set(mv.product_rows()) == ids and rm.product_key(mv) == want:
                return True
        return False

    def channel(self, family, molecules, direction, pathway, product_species) -> RmgChannel:
        keys, orders, statuses = {}, defaultdict(set), Counter()
        directed = 0
        for m in self.moves(family, molecules, direction):
            statuses[m.status] += 1
            if m.status not in ("kept", "unchanged_species") or m.products is None:
                continue
            if m.template is None:
                continue
            inverse = isinstance(pathway, InversePathway)
            if pathway is not None and not inverse and tuple(m.template) != tuple(pathway):
                continue
            if self.species_of_rows(m.product_rows()) != tuple(product_species):
                continue
            if inverse and not self.inverts(family, m, molecules, pathway):
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


class InversePathway(tuple):
    """A pathway filter read through the inverse: the node labels of the
    record's recipe, applied to the moves that undo a counted move (an
    own-reverse family counted on the side opposite to its template)."""


def template_on_products(record) -> bool:
    """Whether the record (or its payload) runs opposite to the reaction RMG
    generated, so its template (the node labels RMG read) sits on its product
    side.  Every count of the record's channel is made on its reactants in its
    direction, whichever member a witness came from, so this, never a
    member's orientation, says whether the count reads the template through
    the inverse (``count_pathway``)."""
    return record.get("orientation", "as_generated") != "as_generated"


def count_pathway(family, pathway, reverse: bool):
    """The pathway filter of the channel count: the record's template.  For a
    family that is its own reverse, counted on the side opposite to its
    template, RMG reads another node on those reactants (the labels are the
    other side's), so the filter is an ``InversePathway``: a counted move
    belongs to the channel when it undoes a move carrying the record's
    template (``RmgCounter.inverts``).  Product species alone would also take
    every other recipe of the family reaching the same species."""
    fam = family
    if getattr(fam, "reverse_template", "x") is None and reverse:
        return InversePathway(pathway) if pathway else None
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
    """(per-temperature max - min of log10 w over the rows, max).  A
    temperature no row measures (no ok row) or measures as non-finite has
    ``None``, and then so has the max: an unmeasured spread is never written
    as inf or NaN, which no payload can carry."""
    per = []
    for i, _t in enumerate(temperatures):
        values = [math.log10(r["w"][i]) for r in rows]
        d = max(values) - min(values) if values else None
        per.append(d if d is not None and math.isfinite(d) else None)
    return per, (max(per) if per and None not in per else None)


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
    reverse: bool = False


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
                               frozenset(r.window.atoms), cert, k == 0, refusal, unextended,
                               reverse=reverse))
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


def candidate_order(c) -> tuple:
    """A total order on candidates read from the candidates themselves: the
    member label, the context index, then the canonical content (witness,
    window, context, certificate, unextended witness, orientation, refusal),
    which also orders members that share a label.  Never the order the
    members come in."""
    def rows(w):
        return None if w is None else [rows_json(w.reactant_rows()), rows_json(w.product_rows())]
    content = [rows(c.witness), sorted(c.window_atoms),
               [[list(p), e] for p, e in c.assignment], c.certificate, rows(c.unextended),
               c.reverse, c.refusal]
    return (c.member, c.context_index, canonical_json(content))


def dedupe(cands) -> list:
    """Candidates with equal executable keys are one template (sources kept).

    An accepted candidate stands for the template whenever the key has one: a
    cut refusal that reads the candidate's own witness (a port reaching a
    physical end) does not refuse an equal executable accepted on another
    witness.  The representative is the least accepted copy in
    ``candidate_order`` (the least refused one when none is accepted), the
    sources are every copy's, sorted, and the templates come base first, then
    by key.  Which templates exist, whether each is refused, which witness
    represents it and the order they are listed in therefore follow the
    candidate set, not the order the members come in (the two records of a
    reciprocal pair draw the same members)."""
    groups, out = {}, []
    for c in cands:
        if c.cut is None:
            c.sources = [(c.member, c.context_index)]
            out.append(c)
        else:
            groups.setdefault(c.cut.key, []).append(c)
    for copies in groups.values():
        rep = min([c for c in copies if c.refusal is None] or copies, key=candidate_order)
        rep.sources = sorted((c.member, c.context_index) for c in copies)
        rep.base = any(c.base for c in copies)
        out.append(rep)
    return sorted(out, key=lambda c: (not c.base, c.cut.key if c.cut else "",
                                      candidate_order(c)))


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
                owner_rate_source=None, kc=None) -> dict:
    """The executable-template payload of one record (``kc``: a ``ThermoKc``;
    when given the payload stores the record's independent thermo Kc).

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
    failing = [{"reason": f"member {case.label[-80:]}: {row.get('status')}"
                          f" {str(row.get('reason', ''))[:160]}",
                "window_status": row.get("status"), "covered": row.get("covered")}
               for case, row, _ctx in members if row.get("status") != "fixed_point"]
    if failing:
        # the member named is the least by its own content, not the first
        return finish({**head, "status": "refused",
                       "refusal": {"by": "upstream", "status": "refused:window",
                                   **min(failing, key=canonical_json)},
                       "templates": [], "contexts": []})
    rec_r = record_graph_rows(record, "reactant_graphs")
    rec_species = counter.species_of_rows(record_graph_rows(record, "product_graphs"))
    rec_reac = counter.species_of_rows(rec_r)
    cands, ctx_rows = [], []
    for case, row, (results, merges) in members:
        if rec_reac == rec_species:
            # a degenerate exchange: one species on both sides cannot say
            # which side the member's witness stands on.  The member was
            # generated with its template side as the case-direction
            # reactants, and the record's template side is its orientation's
            rev = (case.direction != "forward") != template_on_products(record)
        else:
            rev = counter.species_of_rows(dict(case.witness.reactant_rows())) != rec_reac
        on_record_reactants = case.witness.product_rows() if rev else case.witness.reactant_rows()
        if counter.species_of_rows(dict(on_record_reactants)) != rec_reac:
            return finish({**head, "status": "refused",
                           "refusal": {"by": "emission", "status": "refused:orientation",
                                       "reason": "the member's witness carries neither the "
                                                 "record's reactant nor its product species"},
                           "templates": [], "contexts": []})
        c, rows = context_candidates(case, results, merges, family=fam.label, pathway=pathway,
                                     direction=direction, reverse=rev, providers=providers)
        cands.extend(c)
        ctx_rows.extend(dict(r, member=case.label[-120:]) for r in rows)
    cands = dedupe(cands)
    ctx_rows.sort(key=lambda r: (r["member"], r["index"], canonical_json(r)))
    cpath = count_pathway(fam, pathway, template_on_products(record))
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
            # each candidate is counted on its own member's witness, read in
            # that member's orientation (members of both orientations mix),
            # through the record's recipe filter
            cpath=cpath, reverse=cand.reverse,
            counter=counter, rates=rates, weights=weights,
            temperatures=rate["T"], providers=providers, ceiling=ceiling,
            owner_rate_source=owner_rate_source, kc=kc, record_map=next((record_map.get(f"probe{i}") for i, c in enumerate(base)
                             if c is cand), None)))
    status, refusal = record_status(templates)
    out = {**head, "status": status, "refusal": refusal,
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
                   record_map=None, owner_rate_source=None, kc=None) -> dict:
    """Count (item 2, on the context's own witness), locality (item 3) and
    the payload of one candidate."""
    entry = {"member": cand.member[-120:], "context_index": cand.context_index,
             "context": [[p[0], p[1], e] for p, e in cand.assignment if e],
             "base": cand.base, "sources": [[m[-120:], i] for m, i in cand.sources],
             "certificate": cand.certificate,
             # the witnesses below are stored in template direction; this is
             # the side their member was generated on
             "member_reverse": bool(cand.reverse)}
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
    if kc is not None and entry.get("weight"):
        entry["thermo"] = calibration_thermo(kc, counter, rows, weight_source, entry["weight"]["T"])
    # witness graphs carry atom ids: they live in the (unhashed) maps section
    entry["maps"]["witnesses"] = {r["label"]: r.pop("witness") for r in rows if "witness" in r}
    entry["locality"] = {"rows": [record_row] + rows, "spread_log10": per,
                         "max_spread_log10": worst, "tolerance_log10": LOCALITY_TOLERANCE_LOG10,
                         "accepted": ok and worst is not None and
                         worst <= LOCALITY_TOLERANCE_LOG10}
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
    elif worst is None:
        entry.update(status="refused", refusal={
            "status": "refused:locality_spread",
            "reason": "no ok locality row measures the weight" if not good else
                      "a locality row's log10 w is not finite"})
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
               witness={"r": rows_json(r), "p": rows_json(p)},
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
        matcher_from_payload(t, t["weight"]["T"][0], weight=1.0)


def matcher_from_payload(template: dict, temperature: float, *, name=None, primary=None,
                         weight=None) -> wm.Template:
    """The temperature-specific projection of one emitted template to a
    ``window_matcher.Template`` (weight from ``weight_at``: ln w linear in 1/T
    on the grid, a derived reciprocal direction from its ``primary``
    template; extrapolation refused; ``weight`` overrides it when only the
    structure matters)."""
    ex = template["executable"]
    if weight is None:
        weight = weight_at(template, temperature, primary)
    atoms = [wm.WindowAtom(element=a["element"], participant=a["participant"], kind=a["kind"],
                           exterior=tuple(tuple(x) for x in a["exterior"]),
                           isotope=a["isotope"], charge=a["charge"], radical=a["radical"],
                           lone_pairs=a["lone_pairs"]) for a in ex["atoms"]]
    bonds = [(i, j, x) for i, j, x, _y in ex["bonds"] if x]
    edits = [wm.BondEdit(i, j, y) for i, j, x, y in ex["bonds"] if x != y]
    edits += [wm.AtomEdit(i, **a["final"]) for i, a in enumerate(ex["atoms"]) if a["final"]]
    return wm.Template(name or template["name"], template.get("family", "") or "", atoms,
                       bonds, edits, weight=weight)


def weight_at(template: dict, temperature: float, primary=None) -> float:
    """The weight a template loads at ``temperature``.

    Plain templates: ln w linear in 1/T over the grid (``_interpolate``).  A
    template with ``derived`` is the second direction of a reciprocal unit:
    its weight is ``w_primary(T) * (N_P / N_D) / (ssa_ratio * Kc(T))`` with
    ``w_primary`` interpolated as above, ``Kc(T)`` evaluated exactly at ``T``
    from the stored thermo form (``kc_at``), so detailed balance holds at
    every temperature, not only on the grid.  Outside the template's own grid:
    refused."""
    t = float(temperature)
    own = template["weight"]
    d = template.get("derived")
    if d is None:
        return _interpolate(own["T"], own["w"], t)
    if primary is None:
        raise EmissionError("a derived weight needs its primary template")
    ts = [float(x) for x in own["T"]]
    if t < ts[0] or t > ts[-1]:
        raise EmissionError(f"temperature {t} K outside the weight table (extrapolation refused)")
    wp = _interpolate(primary["weight"]["T"], primary["weight"]["w"], t)
    return wp * d["counts_ratio"] / (d["ssa_ratio"] * kc_at(d["thermo"], t))


def _interpolate(ts, ws, t):
    """ln w linear in 1/T between the grid points around ``t``."""
    ts = [float(x) for x in ts]
    if not ts or t < ts[0] or t > ts[-1]:
        raise EmissionError(f"temperature {t} K outside the weight table (extrapolation refused)")
    for i in range(len(ts) - 1):
        if ts[i] <= t <= ts[i + 1]:
            if ts[i + 1] == ts[i]:
                return ws[i]
            f = (1.0 / t - 1.0 / ts[i]) / (1.0 / ts[i + 1] - 1.0 / ts[i])
            return math.exp(math.log(ws[i]) * (1 - f) + math.log(ws[i + 1]) * f)
    return ws[-1]


def interpolation_error(ts, ws) -> tuple:
    """(worst relative error, temperature): each interior grid point predicted
    from its two neighbours by ``_interpolate`` (twice the grid spacing, so an
    upper bound on the error at the grid midpoints) against its own value."""
    worst, at = 0.0, None
    for i in range(1, len(ts) - 1):
        guess = _interpolate([ts[i - 1], ts[i + 1]], [ws[i - 1], ws[i + 1]], float(ts[i]))
        err = abs(guess / ws[i] - 1.0)
        if err > worst:
            worst, at = err, float(ts[i])
    return worst, at


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


def _count_of(payload, template) -> int:
    """``T_ch`` behind a template's weight on its source record's witness."""
    if payload.get("weight_source") == "record":
        n = ((payload.get("counts") or {}).get("record") or {}).get("matcher")
    else:
        n = ((template.get("counts") or {}).get("emitted") or {}).get("matcher")
    if n is None:
        n = ((template.get("counts") or {}).get("context") or {}).get("matcher")
    return int(n) if n else 1


def _pair_deviation(pf, tf, pr, tr, rep_wf, rep_wr) -> tuple:
    """(max over the grid of |log10 (loaded F total / loaded R total) -
    log10 (expected)|, source of the expectation) for one member record pair
    when the loaded (representative) weights replace its own.  Expected: the
    independent thermo ``ssa_ratio * Kc`` of the member's forward template
    (stored by ``emit_record(kc=...)`` on its calibration witness); without it,
    the pair's own calibrated ratio."""
    nf, nr = _count_of(pf, tf), _count_of(pr, tr)
    thermo = tf.get("thermo")
    if thermo is not None and ("error" in thermo or len(thermo["Kc"]) != len(rep_wf)):
        thermo = None
    worst = 0.0
    for i, (lf, lr) in enumerate(zip(rep_wf, rep_wr)):
        loaded = nf * lf / (nr * lr)
        if thermo:
            want = thermo["ssa_ratio"] * thermo["Kc"][i]
        else:
            want = nf * tf["weight"]["w"][i] / (nr * tr["weight"]["w"][i])
        worst = max(worst, abs(math.log10(loaded / want)))
    return worst, "thermo (calibration witness)" if thermo else "own calibrated ratio"


# Where a record id lives (schema slots), masked whole as ``evt_*``: on the
# payload ``record_event_id``, ``reverse_of``, ``rate.provenance.event_id``,
# ``carried_for`` (carry_owners); on a template ``merged_into.record``,
# ``derived.primary.record``, ``collision.records[]``,
# ``collision.pair_deviations[].forward/.reverse``, ``reverse.records[]``,
# ``refusal.records[]`` and ``refusal.weights[].record`` (resolve_collisions),
# ``locality.rows[].provenance.event_id`` (the rate provenance of each
# witness row, ``rate_of``: the id of the record the rate was taken from).
PAYLOAD_ID_SLOTS = (("record_event_id",), ("reverse_of",), ("rate", "provenance", "event_id"),
                    ("carried_for",))
TEMPLATE_ID_SLOTS = (("merged_into", "record"), ("derived", "primary", "record"),
                     ("collision", "records", "*"),
                     ("collision", "pair_deviations", "*", "forward"),
                     ("collision", "pair_deviations", "*", "reverse"),
                     ("reverse", "records", "*"), ("refusal", "records", "*"),
                     ("refusal", "weights", "*", "record"),
                     ("locality", "rows", "*", "provenance", "event_id"))
# Free text that quotes record ids: the payload's and each template's
# ``refusal.reason``.  Only the forms the code writes are masked there: the
# whole id (``evt_`` + 64 hex, e.g. certificates' partner refusals) and its
# 16-character ``str(id)[:16]`` (check_cycles and one_way_problems reasons),
# and only for the ids found in the id slots above.
FREE_TEXT_SLOTS = (("refusal", "reason"),)
ID_TRUNCATIONS = (16,)
# payload fields that are digests over content holding record ids: ``sha256``
# (``finish``: the whole payload, ids included) and ``certificate`` (the
# record certificate's sha256; its owner block holds the record's event id,
# ``certificates.record_owner``); the rate selection it certifies is stored
# in ``rate.provenance.selection``, which the identity keeps
ID_DIGEST_FIELDS = ("sha256", "certificate")
MASK = "evt_*"


def _slot_values(value, path):
    """(container, key) of every place ``path`` names in ``value`` (``*``: each
    list element)."""
    if not path:
        return []
    head, rest = path[0], path[1:]
    if head == "*":
        items = list(enumerate(value)) if isinstance(value, list) else []
    else:
        items = [(head, value[head])] if isinstance(value, dict) and head in value else []
    out = []
    for key, item in items:
        out += [(value, key)] if not rest else _slot_values(item, rest)
    return out


def _mask_ids(value, slots, ids):
    for path in slots:
        for container, key in _slot_values(value, path):
            if isinstance(container[key], str):
                ids.add(container[key])
                container[key] = MASK


def _mask_free_text(value, ids):
    forms = sorted({f for rid in ids if rid != MASK for f in
                    [rid] + [rid[:n] for n in ID_TRUNCATIONS if len(rid) > n]},
                   key=len, reverse=True)
    for path in FREE_TEXT_SLOTS:
        for container, key in _slot_values(value, path):
            text = container[key]
            if isinstance(text, str):
                for form in forms:
                    text = text.replace(form, MASK)
                container[key] = text


def _without_ids(payload, template=None) -> tuple:
    """Copies of ``payload`` (``ID_DIGEST_FIELDS`` dropped) and ``template``
    with the id slots masked whole and, in the free-text slots, the ids found
    there (whole or ``ID_TRUNCATIONS``) masked; every other string as stored."""
    body = json.loads(json.dumps({k: v for k, v in payload.items() if k not in ID_DIGEST_FIELDS}))
    ids = set()
    _mask_ids(body, PAYLOAD_ID_SLOTS, ids)
    for t in body.get("templates", []) if isinstance(body.get("templates"), list) else []:
        _mask_ids(t, TEMPLATE_ID_SLOTS, ids)
    tpl = None
    if template is not None:
        tpl = json.loads(json.dumps(template))
        _mask_ids(tpl, TEMPLATE_ID_SLOTS, ids)
    _mask_free_text(body, ids)
    for t in body.get("templates", []) if isinstance(body.get("templates"), list) else []:
        _mask_free_text(t, ids)
    if tpl is not None:
        _mask_free_text(tpl, ids)
    return body, tpl


def _template_identity(record_identity: str, payload, template) -> str:
    return sha256_json({"payload": record_identity, "template": _without_ids(payload, template)[1]})


def partition_free_identity(payload, template=None) -> str:
    """sha256 of a payload (with ``template``: of the payload and that
    template) as stored, less its record ids and the fields derived from them.

    A record event id hashes the artifact's provenance, the artifact's family
    list among it, so the same chemistry compiled alone or with other
    families carries other ids.  Every choice that must not depend on how the
    families were partitioned (the representative of a collision, the source
    pair of a reciprocal unit, the cycle-check temperatures) orders by this.
    It is the canonical JSON of the whole payload and template with
    ``ID_DIGEST_FIELDS`` dropped and the record ids masked where the schema
    puts them (``_without_ids``), so two copies tie only when their stored
    content is identical, and then any choice between them gives the same
    result."""
    record = sha256_json(_without_ids(payload)[0])
    return record if template is None else _template_identity(record, payload, template)


def resolve_collisions(payloads, tolerance_log10: float = LOCALITY_TOLERANCE_LOG10) -> tuple:
    """One catalog's payloads with equal executable keys made one template.

    Copies on one temperature grid whose log10 weights have a DIAMETER (max -
    min over all copies, at each grid temperature) of at most
    ``tolerance_log10`` can merge, whatever the input order.  The
    representative is the copy of the least ``partition_free_identity`` (then
    the least record event id; never the id alone, which carries the
    artifact's family list); its weight is the
    loaded weight and every source record is listed (``collision``); the
    others become ``merged`` (``merged_into``).  Copies beyond the diameter
    are all refused ``refused:collision``; copies on different grids
    ``refused:collision_grid``.

    Reciprocal pairs (``reverse_of``) are admitted or refused TOGETHER: the
    forward group and its reverse group form one unit that merges only if
    both diameters hold and every member record pair, with the representative
    weights loaded in place of its own, stays within ``tolerance_log10``
    decades of its thermo expectation (``_pair_deviation``; every member's
    deviation is in the table).  The representatives of a unit come from the
    SAME source pair (the pair with the least sorted pair of
    ``partition_free_identity``, then of record ids), so the
    loaded forward and reverse weights are one record pair's.  A partition that
    maps forward copies to more than one reverse template is a
    ``ReciprocalSplitError``.  A non-finite or non-positive weight, or
    ``T_ch = 0``, is an ``EmissionError``.  Returns (new payloads, collision
    table: one row per key with more than one copy)."""
    for p in payloads:
        why = _weight_problem(p)
        if why:
            raise EmissionError(f"record {str(p.get('record_event_id'))[:16]}: {why}")
    out = [json.loads(canonical_json(p)) for p in payloads]
    groups = defaultdict(list)
    identity = {}
    for i, p in enumerate(out):
        record_identity = None
        for j, t in enumerate(p.get("templates", [])):
            if t.get("status") == "emitted":
                groups[t["key"]].append((p.get("record_event_id") or "", i, j))
                record_identity = record_identity or partition_free_identity(p)
                identity[(i, j)] = _template_identity(record_identity, p, t)
    for key in groups:
        groups[key].sort(key=lambda c: (identity[(c[1], c[2])], c))
    links = _reciprocal_partitions(out, groups)
    partner_key = {}
    for (_rid, key), (_rev, rkey) in links.items():
        partner_key.setdefault(key, rkey)
    table, touched, seen = [], set(), set()
    plans = []
    for key in sorted(groups):
        if key in seen:
            continue
        unit = [key]
        other = partner_key.get(key)
        if other is not None and other in groups and other != key:
            unit.append(other)
        seen.update(unit)
        plans.append(unit)
    for unit in plans:
        evals = {}
        for key in unit:
            copies = groups[key]
            weights = [out[i]["templates"][j]["weight"] for _r, i, j in copies]
            grid = next((_grid_problem(weights[0]["T"], w["T"]) for w in weights[1:]
                         if _grid_problem(weights[0]["T"], w["T"])), None)
            diameter = None
            if grid is None:
                logs = [[math.log10(x) for x in w["w"]] for w in weights]
                diameter = max(max(c) - min(c) for c in zip(*logs))
            evals[key] = {"copies": copies, "weights": weights, "grid": grid,
                          "diameter": diameter}
        if all(len(evals[k]["copies"]) < 2 for k in unit):
            continue
        reps = {k: evals[k]["copies"][0] for k in unit}
        devs, dev_rows = 0.0, []
        if len(unit) == 2:
            kf, kr = unit
            ident = {r: identity[(i, j)] for k in unit for r, i, j in evals[k]["copies"]}
            pairs = sorted(((f_rid, links[(f_rid, kf)][0]) for f_rid, _i, _j in evals[kf]["copies"]
                            if (f_rid, kf) in links),
                           key=lambda x: (tuple(sorted(ident.get(r, "") for r in x)), tuple(sorted(x))))
            if pairs:
                rep_pair = pairs[0]
                reps[kf] = next(c for c in evals[kf]["copies"] if c[0] == rep_pair[0])
                reps[kr] = next(c for c in evals[kr]["copies"] if c[0] == rep_pair[1])
                if all(e["grid"] is None for e in evals.values()):
                    rf = out[reps[kf][1]]["templates"][reps[kf][2]]["weight"]["w"]
                    rr = out[reps[kr][1]]["templates"][reps[kr][2]]["weight"]["w"]
                    by_rid = {r: (i, j) for k in unit for r, i, j in evals[k]["copies"]}
                    for f_rid, r_rid in pairs:
                        fi, fj = by_rid[f_rid]
                        ri, rj = by_rid[r_rid]
                        dev, how = _pair_deviation(
                            out[fi], out[fi]["templates"][fj], out[ri], out[ri]["templates"][rj],
                            rf, rr)
                        dev_rows.append({"forward": f_rid, "reverse": r_rid,
                                         "deviation_log10": dev, "expected": how})
                        devs = max(devs, dev)
        grid_bad = next((k for k in unit if evals[k]["grid"]), None)
        diam_bad = [k for k in unit if evals[k]["diameter"] is not None
                    and beyond_tolerance(evals[k]["diameter"], tolerance_log10)]
        if grid_bad is not None:
            admit, disposition = False, COLLISION_STATUS + "_grid"
            reason = f"temperature grids differ: {evals[grid_bad]['grid']}"
        elif diam_bad or beyond_tolerance(devs, tolerance_log10):
            admit, disposition = False, COLLISION_STATUS
            reason = (f"log10 weight diameter {[evals[k]['diameter'] for k in unit]} or loaded "
                      f"pair deviation {devs:.4g} > {tolerance_log10}")
        else:
            admit, disposition, reason = True, "merged", ""
        for key in unit:
            ev = evals[key]
            copies = ev["copies"]
            records = [r for r, _i, _j in copies]
            spread = None
            if ev["grid"] is None:
                spread = max((max(c) - min(c)) / max(c)
                             for c in zip(*[w["w"] for w in ev["weights"]]))
            if len(copies) >= 2 or len(unit) == 2 and not admit:
                table.append({"key": key, "records": records, "copies": len(copies),
                              "log10_diameter": ev["diameter"], "rel_spread": spread,
                              "disposition": disposition, "tolerance_log10": tolerance_log10,
                              "unit": [k[:16] for k in unit],
                              "pair_deviations": dev_rows if len(unit) == 2 else [],
                              "max_pair_deviation_log10": devs if len(unit) == 2 else None})
            if len(copies) < 2:
                continue
            owner_record, oi, oj = reps[key]
            weights = ev["weights"]
            for r, i, j in copies:
                t = out[i]["templates"][j]
                touched.add(i)
                if admit and (i, j) == (oi, oj):
                    t["collision"] = {
                        "merged": True, "records": records,
                        "log10_diameter": ev["diameter"], "rel_spread": spread,
                        "tolerance_log10": tolerance_log10,
                        "weight_rule": "the representative's weight: the copy (of a unit: the "
                                       "source pair) of the least partition-free identity, "
                                       "record ids breaking ties only",
                        "pair_deviations": dev_rows if len(unit) == 2 else []}
                elif admit:
                    t["status"] = "merged"
                    t["merged_into"] = {"record": owner_record, "key": key}
                else:
                    t["status"] = "refused"
                    t["refusal"] = {"status": disposition, "reason": reason[:300], "keys": [key],
                                    "records": records,
                                    "weights": [{"record": rr, "T": w["T"], "w": w["w"]}
                                                for rr, w in zip(records, weights)]}
    live = defaultdict(set)
    for p in out:
        for t in p.get("templates", []):
            if t.get("status") in ("emitted", "merged") and t.get("base"):
                rid = p.get("record_event_id")
                live[t["merged_into"]["key"] if t["status"] == "merged" else t["key"]].add(rid)
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
    # reciprocal groups were decided together: a unit half merged is a defect
    final = {}
    for p in out:
        for t in p.get("templates", []):
            if not t.get("base") or t["status"] == "refused" and "key" not in t:
                continue          # refused at emission: never cut, so in no unit
            final[t["key"]] = t["status"] if t["status"] != "merged" else "emitted"
    for (_rid, key), (_rev, rkey) in links.items():
        a, b = final.get(key), final.get(rkey)
        if a is not None and b is not None and (a == "refused") != (b == "refused"):
            raise ReciprocalSplitError(
                f"reciprocal split: template {key[:16]} is {a} but its reverse {rkey[:16]} is {b}")
    for i in touched:
        out[i]["status"], out[i]["refusal"] = record_status(out[i]["templates"])
        out[i] = finish(out[i])
    return out, table


# ---------------------------------------------------------------------------
# independent thermo and the equilibrium cycle on the loaded catalog
# ---------------------------------------------------------------------------
def rows_json(rows) -> list:
    return [[a, x[0], x[1], x[2], x[3], x[4], [[n, o] for n, o in sorted(x[5])]]
            for a, x in sorted(dict(rows).items())]


def rows_from_json(items) -> dict:
    return {a: (el, iso, ch, rad, lp, tuple((n, float(o)) for n, o in edges))
            for a, el, iso, ch, rad, lp, edges in items}


def ssa_of_rows(rows, counter) -> float:
    """The pair convention of a reactant side: 2 for two copies of one species
    (the SSA counts N(N-1)/2 pairs), else 1 (``event_record.ssa_multiplier_for``)."""
    parts = rm.split_rows(dict(rows))
    if len(parts) == 2 and counter.species_of_rows(parts[0]) == counter.species_of_rows(parts[1]):
        return 2.0
    return 1.0


KC_FORM = "ThermoData-Kc/1"
_THERMO_CACHE = {}


def _thermo_json(thermo) -> dict:
    """The numbers of an RMG ``ThermoData`` (the form every species of the
    reference thermo path takes); anything else is not storable."""
    from rmgpy.thermo import ThermoData
    if not isinstance(thermo, ThermoData):
        raise EmissionError(f"thermo of type {type(thermo).__name__} cannot be stored "
                            f"(only ThermoData)")
    return {"Tdata": [float(x) for x in thermo.Tdata.value_si],
            "Cpdata": [float(x) for x in thermo.Cpdata.value_si],
            "H298": float(thermo.H298.value_si), "S298": float(thermo.S298.value_si),
            "Cp0": float(thermo.Cp0.value_si), "CpInf": float(thermo.CpInf.value_si)}


def _thermo_object(item):
    key = canonical_json(item)
    if key not in _THERMO_CACHE:
        from rmgpy.thermo import ThermoData
        _THERMO_CACHE[key] = ThermoData(
            Tdata=(item["Tdata"], "K"), Cpdata=(item["Cpdata"], "J/(mol*K)"),
            H298=(item["H298"], "J/mol"), S298=(item["S298"], "J/(mol*K)"),
            Cp0=(item["Cp0"], "J/(mol*K)"), CpInf=(item["CpInf"], "J/(mol*K)"))
    return _THERMO_CACHE[key]


def kc_at(form: dict, temperature: float) -> float:
    """``Kc(T)`` of a stored thermo form, exactly at ``temperature``: the
    formula of ``Reaction.get_equilibrium_constant(T, type="Kc")`` (``Ka =
    exp(-dG/RT)``, ``Kc = Ka * C0**(n_products - n_reactants)``, ``C0 = P0/RT``,
    ``P0 = 1e5 Pa``) with each species' free energy from the stored
    ``ThermoData``."""
    from rmgpy import constants
    if form.get("kind") != KC_FORM:
        raise EmissionError(f"unknown thermo form {form.get('kind')!r}")
    t = float(temperature)
    g = lambda side: math.fsum(_thermo_object(x).get_free_energy(t) for x in side)
    dg = g(form["products"]) - g(form["reactants"])
    c0 = form.get("P0", 1e5) / constants.R / t
    return float(math.exp(-dg / constants.R / t) * c0 ** (len(form["products"])
                                                         - len(form["reactants"])))


class ThermoKc:
    """Whole-species ``Kc(T)`` (concentration units of the displayed reaction
    ``reactants -> products``) from RMG gas-phase thermo, the path the compile
    used: ``reference_thermo.GasPhaseRMGReferenceThermo.evaluate`` on a
    ``Reaction`` of one ``Species`` per connected component, thermo assigned by
    ``SharedThermoAssignment`` (compiler.py:4268-4277 builds the same provider).

    ``form(r, p, T)`` is the stored representation (``KC_FORM``: the
    ``ThermoData`` numbers of every participant) that ``kc_at`` evaluates at
    any temperature.

    Detailed-balance convention (see ``check_cycles``): the loaded totals of
    the two directions satisfy ``S_f / S_r = (ssa_f / ssa_r) * Kc``, with
    ``ssa`` 2 for an identical-species reactant pair (``ssa_of_rows``)."""

    def __init__(self, provider):
        self.provider = provider
        self._cache = {}
        self._forms = {}

    def _reaction(self, reactant_rows, product_rows):
        from rmgpy.reaction import Reaction
        from rmgpy.species import Species

        def side(rows):
            return [Species(molecule=[rm.molecule_from_rows(part)])
                    for part in rm.split_rows(dict(rows))]
        return Reaction(reactants=side(reactant_rows), products=side(product_rows))

    def __call__(self, reactant_rows, product_rows, temperatures) -> list:
        key = (rm._signature(dict(reactant_rows)), rm._signature(dict(product_rows)),
               tuple(temperatures))
        if key not in self._cache:
            reaction = self._reaction(reactant_rows, product_rows)
            try:
                values = self.provider.evaluate(reaction, list(temperatures)).equilibrium_constants
            except ValueError as error:      # ThermoUnavailable
                raise EmissionError(f"thermo Kc unavailable: {error}") from None
            self._cache[key] = [float(v) for v in values]
            try:
                self._forms[key[:2]] = {
                    "kind": KC_FORM, "P0": 1e5,
                    "reactants": [_thermo_json(x.thermo) for x in reaction.reactants],
                    "products": [_thermo_json(x.thermo) for x in reaction.products]}
            except EmissionError as error:
                self._forms[key[:2]] = {"kind": "unstorable", "reason": str(error)}
        return list(self._cache[key])

    def form(self, reactant_rows, product_rows, temperatures) -> dict:
        self(reactant_rows, product_rows, temperatures)
        return self._forms[(rm._signature(dict(reactant_rows)),
                            rm._signature(dict(product_rows)))]


def calibration_thermo(kc, counter, rows, weight_source, grid) -> dict:
    """The thermo of the witness a template's weight was calibrated on (the
    record's own witness, or the emitted witness for a weight that reaches
    beyond the proxy graph): the stored form, its grid Kc and the ``ssa``
    ratio, in the template direction."""
    label = "context" if weight_source == "record" else "emitted"
    row = next((r for r in rows if r["label"] == label and r.get("status") == "ok"
                and "witness" in r), None)
    if row is None:
        return {"error": f"no {label} witness to take the thermo from"}
    r_rows, p_rows = rows_from_json(row["witness"]["r"]), rows_from_json(row["witness"]["p"])
    form = kc.form(r_rows, p_rows, grid)
    if form.get("kind") != KC_FORM:
        return {"error": form.get("reason", "unstorable thermo")}
    return {"witness": label, "T": list(grid), "Kc": list(kc(r_rows, p_rows, grid)),
            "ssa_ratio": ssa_of_rows(r_rows, counter) / ssa_of_rows(p_rows, counter),
            "form": form,
            "convention": "S_f/S_r = ssa_ratio * Kc (template direction, calibration witness)"}


class _Loaded:
    """The templates a catalog loads (``emitted``), by (family, direction),
    with the weights they have at any requested temperature."""

    def __init__(self, payloads):
        self.by = defaultdict(list)
        self.index = {}
        self.owner = {}
        for p in payloads:
            for t in p.get("templates", []):
                if t.get("status") == "emitted":
                    self.by[(t.get("family") or p.get("family"), t["direction"])].append(t)
                    self.index[(p.get("record_event_id"), t["key"])] = t
                    self.owner[id(t)] = (p.get("record_event_id"), t["key"])

    def weight(self, template, temperature):
        d = template.get("derived")
        primary = None
        if d is not None:
            primary = self.index.get((d["primary"]["record"], d["primary"]["key"]))
            if primary is None:
                raise EmissionError("the primary template of a derived weight is not loaded")
        return weight_at(template, temperature, primary)

    def total(self, counter, rows, payload, direction, target_rows, temperatures):
        """(weights at ``temperatures`` summed over the loaded transitions of
        ``payload``'s channel on ``rows``, the transitions, the state)."""
        fam = counter.family(payload["family"])
        channel = counter.channel(
            fam, molecules_of(rows), direction,
            count_pathway(fam, tuple(payload.get("pathway") or ()),
                          template_on_products(payload)),
            counter.species_of_rows(target_rows))
        keys = set(channel.keys)
        entries = self.by.get((fam.label if hasattr(fam, "label") else payload["family"],
                               direction), [])
        templates = [matcher_from_payload(dict(t, family=payload["family"]), 0.0, weight=1.0,
                                          name=f"L{k}") for k, t in enumerate(entries)]
        state = melt_state(rows)
        inside = []
        if templates:
            inside = [tr for tr in wm.match(templates, state).transitions
                      if rm.product_key(wm.product_rows(state, tr)) in keys]
        sums = [math.fsum(self.weight(entries[int(tr.template[1:])], T) for tr in inside)
                for T in temperatures]
        return sums, inside, state


def _norm_row(row):
    return (row[0], row[1], int(row[2]), int(row[3]), int(row[4]),
            tuple(sorted((int(n), float(o)) for n, o in row[5])))


def _inverse_counts(loaded, counter, rows, state, inside, payload, direction):
    """Per forward transition: how many loaded transitions of ``payload`` on
    its product state give the original components back (1 expected)."""
    entries = loaded.by.get((payload["family"], direction), [])
    templates = [matcher_from_payload(dict(t, family=payload["family"]), 0.0, weight=1.0,
                                      name=f"L{k}") for k, t in enumerate(entries)]
    first = first_form_rows(rows)
    counts = []
    for tr in inside:
        touched = wm.product_rows(state, tr)
        after = dict(rows)
        after.update(touched)
        check = melt_state(after)
        back = 0
        for rt in wm.match(templates, check).transitions if templates else []:
            got = wm.product_rows(check, rt)
            if set(got) == set(touched) and \
                    all(_norm_row(got[a]) == _norm_row(first[a]) for a in got):
                back += 1
        counts.append(back)
    return counts


CHECK_LABELS = ("context", "emitted")


BOUNDARY_SLACK_LOG10 = 1e-12      # float slack: 10**0.01 scaling measures 0.01 + 6.8e-16


def beyond_tolerance(distance_log10: float, tolerance_log10: float = LOCALITY_TOLERANCE_LOG10) -> bool:
    """The one boundary of "the same weight": two weights (or a weight and its
    expectation) are at ``|log10 ratio|`` decades; beyond ``tolerance_log10``
    (default ``LOCALITY_TOLERANCE_LOG10``, the collision diameter's) is
    refused, at it (within ``BOUNDARY_SLACK_LOG10`` of float noise) accepted.  The merged-copy guard, the merged
    pair's ratio check and the collision diameter all use this comparison, so
    no ratio is accepted by one and rejected by another (a linear |ratio - 1|
    bound is asymmetric in the ratio)."""
    return distance_log10 > tolerance_log10 + BOUNDARY_SLACK_LOG10


def _loads(payload) -> bool:
    """The record loads a template (its own, or its representative's)."""
    return any(t.get("status") in ("emitted", "merged") for t in payload.get("templates", []))


def one_way_problems(payloads, complete: bool = True) -> list:
    """(record, reason) for every record that loads a template while its
    ``reverse_of`` record is absent (``complete`` catalogs) or loads none (a
    one-way channel)."""
    by_record = {p.get("record_event_id"): p for p in payloads}
    out = []
    for rid, p in sorted(by_record.items(), key=lambda x: str(x[0])):
        rev = p.get("reverse_of")
        if rev is None or not _loads(p):
            continue
        q = by_record.get(rev)
        if q is None:
            if complete:
                out.append((rid, f"reverse record {str(rev)[:16]} is not in the catalog"))
        elif not _loads(q):
            out.append((rid, f"reverse record {str(rev)[:16]} loads no template"))
    return out


def _inverse_partners(loaded, counter, rows, state, inside, payload, direction) -> set:
    """The loaded templates (record, key) of ``payload``'s family and
    ``direction`` whose transition on a forward transition's product state
    gives the original components back."""
    entries = loaded.by.get((payload["family"], direction), [])
    templates = [matcher_from_payload(dict(t, family=payload["family"]), 0.0, weight=1.0,
                                      name=f"L{k}") for k, t in enumerate(entries)]
    first = first_form_rows(rows)
    found = set()
    for tr in inside:
        touched = wm.product_rows(state, tr)
        after = dict(rows)
        after.update(touched)
        check = melt_state(after)
        for rt in wm.match(templates, check).transitions if templates else []:
            got = wm.product_rows(check, rt)
            if set(got) == set(touched) and \
                    all(_norm_row(got[a]) == _norm_row(first[a]) for a in got):
                found.add(loaded.owner[id(entries[int(rt.template[1:])])])
    return found


def primary_order(payload, template, record_id) -> tuple:
    """The order of the two templates of a reciprocal pair; the lesser is the
    primary (``derive_weights``): as_generated before reversed, then the
    template key, then (never decisive between loaded templates, whose keys
    are unique) the record event id."""
    return (payload.get("orientation") != "as_generated", template.get("key"), str(record_id))


def derive_weights(payloads, *, counter) -> list:
    """Make the weights of every reciprocal pair of LOADED templates follow
    detailed balance at the requested temperature.

    For each emitted template ``t`` (base or context) with a stored thermo, the
    reverse templates that undo its transitions on its ``context`` and
    ``emitted`` witnesses are found by matching the loaded catalog.  When
    there is exactly one partner ``u``, the pair is ordered by
    ``primary_order``: the template of the record RMG generated (orientation
    ``as_generated``, the record whose own k the compiler keeps, deriving its
    partner by Kc) first, then the lesser template key (the sha256 of the
    executable).  Both are the same however the families were partitioned; a
    record id hashes the artifact's family list and only breaks a tie, which
    loaded keys never have.  The lesser is PRIMARY and keeps its interpolated
    table, the greater is DERIVED: ``derived`` names the primary and carries the
    primary's stored thermo form, ``N_P / N_D`` and the ``ssa`` ratio, so that
    ``weight_at`` returns ``w_P(T) * (N_P / N_D) / (ssa_ratio * Kc_P(T))``.  A
    template with several or no partners, or without thermo, stays
    independent (and fails the off-grid cycle check)."""
    out = [json.loads(canonical_json(p)) for p in payloads]
    loaded = _Loaded(out)
    by_record = {p.get("record_event_id"): p for p in out}
    changed = set()
    for i, p in enumerate(out):
        rid, rev = p.get("record_event_id"), p.get("reverse_of")
        q = by_record.get(rev)
        if q is None:
            continue
        for t in p.get("templates", []):
            if t.get("status") != "emitted" or t.get("derived"):
                continue
            partners = set()
            for lab in CHECK_LABELS:
                wit = ((t.get("maps") or {}).get("witnesses") or {}).get(lab)
                row = next((r for r in (t.get("locality") or {}).get("rows", [])
                            if r.get("label") == lab and r.get("status") == "ok"), None)
                if wit is None or row is None:
                    continue
                r_rows, p_rows = rows_from_json(wit["r"]), rows_from_json(wit["p"])
                try:
                    _s, inside, state = loaded.total(counter, r_rows, p, t["direction"], p_rows, [])
                    entries = loaded.by.get((p["family"], t["direction"]), [])
                    mine = [tr for tr in inside if entries[int(tr.template[1:])] is t]
                    partners |= _inverse_partners(loaded, counter, r_rows, state, mine, q,
                                                  q["direction"])
                except EmissionError:
                    partners.add(None)
            if len(partners) != 1 or None in partners:
                continue
            ref = next(iter(partners))
            if ref == (rid, t["key"]):
                # a self-inverse class (after merging, its own reverse): primary,
                # never derived; its balance is checked directly
                t["self_inverse"] = True
                changed.add(i)
                continue
            u = loaded.index[ref]
            if primary_order(p, t, rid) < primary_order(by_record.get(ref[0]) or {}, u, ref[0]):
                continue                         # t is the primary of this pair
            th = u.get("thermo")
            if not th or "error" in th or "form" not in th:
                continue
            t["derived"] = {
                "primary": {"record": ref[0], "key": ref[1]},
                "counts_ratio": _count_of(by_record[ref[0]], u) / _count_of(p, t),
                "ssa_ratio": th["ssa_ratio"], "thermo": th["form"],
                "rule": "primary = of the pair of loaded templates that undo each other, the "
                        "one of the as_generated record, then of the lesser template key (record "
                        "ids, which carry the artifact's family list, only break ties); "
                        "w = w_primary(T) * (N_P / N_D) / (ssa_ratio * Kc(T)), Kc exactly at T"}
            changed.add(i)
    return [finish(p) if k in changed else p for k, p in enumerate(out)]


OFFGRID_RANDOM = 5
DERIVED_RTOL = 1e-9


def check_temperatures(grid, pair, extra=()) -> list:
    """Where a loaded pair is checked: every grid temperature, every grid
    midpoint, ``OFFGRID_RANDOM`` in-range temperatures drawn from a generator
    seeded by ``pair`` (``check_cycles`` passes the pair's sorted
    ``partition_free_identity``, never its record ids), and ``extra``."""
    import random
    grid = [float(x) for x in grid]
    mids = [(a + b) / 2.0 for a, b in zip(grid, grid[1:])]
    rng = random.Random(int(hashlib.sha256("|".join(pair).encode()).hexdigest()[:16], 16))
    rand = [grid[0] + (grid[-1] - grid[0]) * rng.random() for _ in range(OFFGRID_RANDOM)]
    return sorted(set(grid + mids + rand + [float(x) for x in extra]))


def check_cycles(payloads, *, counter, kc, complete: bool = True, rtol: float = CYCLE_RTOL,
                 extra_temperatures=(), max_rounds: int = 6) -> tuple:
    """Detailed balance of the LOADED catalog against independent thermo, at
    every temperature the runtime can request.

    Applies to every payload with a ``reverse_of`` record that loads (or
    merges) a template.  For each ``context`` and ``emitted`` locality row of
    each such template, ``emitted`` or ``merged`` (the row stores its
    witness's template-direction reactants ``r`` and products ``p``), at the
    grid temperatures, the grid midpoints, ``OFFGRID_RANDOM`` seeded random
    in-range temperatures per pair and ``extra_temperatures``:

    * ``S_f``: the sum of the weights the LOADED templates have at ``T``
      (``weight_at``; a merged copy loads its representative) over the
      transitions of the record's channel on ``r``;
    * ``S_r``: the same for the reverse record's channel on ``p``;
    * bijection: each forward transition is undone by exactly one loaded
      reverse transition on its product state;
    * ratio: ``S_f / S_r = (ssa(r) / ssa(p)) * Kc(r -> p; T)`` with ``Kc`` from
      ``kc`` (``ThermoKc``: whole-species RMG thermo of the witness species,
      concentration units, never the records' k_F / k_R) and ``ssa`` 2 for an
      identical-species pair.  Tolerance ``rtol`` (1e-12) on an ``emitted``
      template's rows; a ``merged`` copy's rows may deviate by the locality
      tolerance, and the deviation is reported;
    * a derived template's stored (record) weights equal its derivation at
      the grid to ``DERIVED_RTOL`` (a pair whose own kinetics are not
      reverse-by-Kc is refused, not silently replaced).

    Un-evaluable (reverse record absent in a ``complete`` catalog, a missing
    witness or thermo, no loaded reverse transition) is a refusal, not a
    skip; with ``complete=False`` an absent reverse is ``not_evaluable``.  A
    failing pair refuses every emitted template of both records, the
    representatives behind their merged copies, and every merged copy of
    those; a record whose reverse record loads nothing is refused the same
    way; repeated until stable.  Returns (new payloads, one row per
    reciprocal record pair)."""
    out = [json.loads(canonical_json(p)) for p in payloads]
    by_record = {p.get("record_event_id"): i for i, p in enumerate(out)}
    # the seed of each pair's random check temperatures: the records' content
    # as given (fixed over the rounds), never their ids
    seed_of = {rid: partition_free_identity(out[i]) for rid, i in by_record.items()}
    verdict = {}            # pair -> table row (refused ones stick)
    for _round in range(max_rounds):
        loaded = _Loaded(out)
        rows_out, refuse = {}, {}
        cache = {}
        for rid in sorted(by_record, key=str):
            p = out[by_record[rid]]
            rev = p.get("reverse_of")
            if rev is None:
                continue
            pair = tuple(sorted((rid, rev)))
            mine = [t for t in p.get("templates", []) if t.get("status") in ("emitted", "merged")]
            if not mine:
                continue
            row = rows_out.setdefault(pair, {
                "forward": pair[0], "reverse": pair[1], "rows": 0, "max_ratio_error": 0.0,
                "max_grid_error": 0.0, "max_offgrid_error": 0.0, "temperatures": 0,
                "bijection_failures": 0, "status": "ok", "reason": "", "problems": []})
            q = out[by_record[rev]] if rev in by_record else None

            def fail(why, row=row, pair=pair):
                row["status"] = "refused:cycle"
                if why not in row["problems"]:
                    row["problems"].append(why)
                row["reason"] = "; ".join(row["problems"])[:500]
                refuse[pair] = row["reason"]
            if q is None:
                if complete:
                    fail(f"reverse record {str(rev)[:16]} is not in the catalog")
                else:
                    row["status"] = "not_evaluable"
                    row["reason"] = "reverse record outside this sample"
                continue
            if not _loads(q):
                fail(f"reverse record {str(rev)[:16]} loads no template (one-way channel)")
                continue
            seed_pair = tuple(sorted((seed_of[rid], seed_of[rev])))
            for t in mine:
                if t.get("status") == "merged":
                    # a merged copy's own record kinetics must match the loaded
                    # weight (representative's, or derived) within the locality
                    # tolerance at every grid temperature
                    into = t["merged_into"]
                    owner = loaded.index.get((into["record"], into["key"]))
                    try:
                        if owner is None:
                            raise EmissionError("its owner is not loaded")
                        bad = max(abs(math.log10(loaded.weight(owner, T) / w))
                                  for T, w in zip(t["weight"]["T"], t["weight"]["w"]))
                        if beyond_tolerance(bad):
                            fail(f"merged template {t['key'][:12]}: the record's own kinetics "
                                 f"differ from the loaded weight by {bad:.3g} decades > "
                                 f"{LOCALITY_TOLERANCE_LOG10:.3g}")
                    except EmissionError as error:
                        fail(f"merged template {t['key'][:12]}: {error}")
                if t.get("status") == "emitted" and t.get("derived"):
                    try:
                        grid = t["weight"]["T"]
                        got = [loaded.weight(t, T) for T in grid]
                        bad = max(abs(a / b - 1.0) for a, b in zip(got, t["weight"]["w"]))
                        if bad > DERIVED_RTOL:
                            fail(f"template {t['key'][:12]}: the record's own kinetics differ "
                                 f"from the thermo-derived weight by {bad:.3g} (not a "
                                 f"reverse-by-Kc pair)")
                    except EmissionError as error:
                        fail(f"template {t['key'][:12]}: {error}")
                for lr in (t.get("locality") or {}).get("rows", []):
                    if lr.get("label") not in CHECK_LABELS or lr.get("status") != "ok":
                        continue
                    wit = ((t.get("maps") or {}).get("witnesses") or {}).get(lr["label"])
                    if wit is None:
                        fail(f"template {t['key'][:12]} row {lr['label']}: no stored witness")
                        continue
                    r_rows = rows_from_json(wit["r"])
                    p_rows = rows_from_json(wit["p"])
                    temps = check_temperatures(t["weight"]["T"], seed_pair, extra_temperatures)
                    ck = (rm._signature(r_rows), rm._signature(p_rows), t["direction"],
                          p["family"], q["family"], tuple(p.get("pathway") or ()),
                          tuple(q.get("pathway") or ()), tuple(temps))
                    if ck not in cache:
                        try:
                            s_f, inside, state = loaded.total(
                                counter, r_rows, p, t["direction"], p_rows, temps)
                            s_r, inside_r, _st = loaded.total(
                                counter, p_rows, q, q["direction"], r_rows, temps)
                            if not inside or not inside_r:
                                raise EmissionError(
                                    f"no loaded {'forward' if not inside else 'reverse'} "
                                    f"transition on the witness")
                            back = _inverse_counts(loaded, counter, r_rows, state, inside, q,
                                                   q["direction"])
                            kcs = kc(r_rows, p_rows, temps)
                            want = [ssa_of_rows(r_rows, counter) / ssa_of_rows(p_rows, counter) * k
                                    for k in kcs]
                            errs = [abs(a / b / c - 1.0) for a, b, c in zip(s_f, s_r, want)]
                            decades = [abs(math.log10(a / b / c))
                                       for a, b, c in zip(s_f, s_r, want)]
                            cache[ck] = {"errs": errs, "decades": decades, "back": back, "grid": temps, "s_f": s_f,
                                         "s_r": s_r, "kc": kcs, "want": want}
                        except EmissionError as error:
                            cache[ck] = {"error": str(error)}
                    res = cache[ck]
                    row["rows"] += 1
                    if "error" in res:
                        fail(f"template {t['key'][:12]} {lr['label']}: {res['error']}"[:300])
                        continue
                    on_grid = set(float(x) for x in t["weight"]["T"])
                    row["temperatures"] = max(row["temperatures"], len(res["grid"]))
                    for T, e in zip(res["grid"], res["errs"]):
                        key_ = "max_grid_error" if T in on_grid else "max_offgrid_error"
                        row[key_] = max(row[key_], e)
                    emitted_t = t.get("status") == "emitted"
                    worst = max(res["errs"])
                    i = res["errs"].index(worst)
                    if worst > row["max_ratio_error"] or "worst" not in row:
                        row["max_ratio_error"] = worst
                        row["worst"] = {"T": res["grid"][i], "S_f": res["s_f"][i],
                                        "S_r": res["s_r"][i], "Kc": res["kc"][i],
                                        "expected_ratio": res["want"][i],
                                        "ongrid": res["grid"][i] in on_grid,
                                        "template": t["key"][:16], "status": t["status"],
                                        "row": lr["label"]}
                    if extra_temperatures:
                        row.setdefault("extra", {}).update(
                            {T: e for T, e in zip(res["grid"], res["errs"])
                             if T in set(float(x) for x in extra_temperatures)})
                    if any(n != 1 for n in res["back"]):
                        row["bijection_failures"] += 1
                        fail(f"inverse edit: each forward transition must be undone by exactly "
                             f"one loaded reverse transition (counts {res['back']})")
                    if (worst > rtol) if emitted_t else beyond_tolerance(max(res["decades"])):
                        fail(f"ratio: S_f/S_r = {res['s_f'][i] / res['s_r'][i]:.17g} vs "
                             f"(ssa ratio * Kc) = {res['want'][i]:.17g} at {res['grid'][i]} K"
                             f"{'' if res['grid'][i] in on_grid else ' (off the grid)'}, "
                             f"relative error {worst:.3g} ("
                             + (f"> {rtol:.3g}" if emitted_t else f"{max(res['decades']):.3g} "
                                f"decades > {LOCALITY_TOLERANCE_LOG10:.3g}")
                             + f") (template "
                             f"{t['key'][:12]}, {t['status']}, row {lr['label']})")
        for pair, why in refuse.items():
            verdict[pair] = rows_out[pair]
        for pair, row in rows_out.items():
            if pair not in verdict and row["status"] != "ok":
                verdict[pair] = row
        dead = {}

        def kill(rid, why):
            for t in out[by_record[rid]].get("templates", []):
                if t.get("status") == "emitted":
                    dead[(rid, t["key"])] = why
                elif t.get("status") == "merged":
                    dead[(t["merged_into"]["record"], t["merged_into"]["key"])] = why
        for pair, why in refuse.items():
            for rid in pair:
                if rid in by_record:
                    kill(rid, why)
        for rid, why in one_way_problems(out, complete):
            kill(rid, why)
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
                t["refusal"] = {"status": "refused:cycle", "reason": str(why)[:600]}
                t.pop("merged_into", None)
                changed.add(i)
        for i in sorted(changed):
            out[i]["status"], out[i]["refusal"] = record_status(out[i]["templates"])
            out[i] = finish(out[i])
        if not changed:
            final = {pair: rows_out[pair] for pair in rows_out}
            for pair, row in verdict.items():
                final[pair] = row
            return out, [final[k] for k in sorted(final)]
    raise EmissionError(f"the cycle check did not settle in {max_rounds} rounds")


def cycle_gates(table, payloads, complete: bool = True) -> list:
    """Acceptance gates: every refused pair and every loaded template whose
    reverse record loads nothing.  Pairs a partial catalog cannot evaluate
    (``not_evaluable``) are reported separately (``cycle_not_evaluable``)."""
    gates = [f"cycle {r['forward'][:16]}/{r['reverse'][:16]}: {r['status']}: {r['reason']}"[:500]
             for r in table if r["status"] not in ("ok", "not_evaluable")]
    gates += [f"one-way channel {str(rid)[:16]}: {why}" for rid, why in one_way_problems(payloads, complete)]
    return gates


def cycle_not_evaluable(table) -> list:
    """The pairs a partial catalog could not evaluate (not gates)."""
    return [f"{r['forward'][:16]}/{r['reverse'][:16]}: {r['reason']}"
            for r in table if r["status"] == "not_evaluable"]


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
    index = {(p.get("record_event_id"), t["key"]): t for p in payloads
             for t in p.get("templates", []) if t.get("status") == "emitted"}
    out = []
    for p in payloads:
        for t in p.get("templates", []):
            if t.get("status") == "emitted":
                label = (name or (lambda p, t: f"{(p.get('record_event_id') or '')[4:16]}:"
                                               f"{t['key'][:12]}"))(p, t)
                d = t.get("derived")
                primary = None
                if d is not None:
                    primary = index.get((d["primary"]["record"], d["primary"]["key"]))
                    if primary is None:
                        raise EmissionError(
                            f"orphaned derived template {t['key'][:16]}: its primary "
                            f"{d['primary']['record'][:16]} is not loaded")
                out.append(matcher_from_payload(t, temperature, name=label, primary=primary))
    return out


def carry_owners(selected, catalog) -> list:
    """``selected`` payloads (a subset of ``catalog``) plus, for each merged
    template in them whose owner is not selected, and for each derived
    template whose primary is not selected, a payload of that one template
    (not the other templates of its record)."""
    have = {(p.get("record_event_id"), t.get("key")) for p in selected
            for t in p.get("templates", []) if t.get("status") == "emitted"}
    by_record = {p.get("record_event_id"): p for p in catalog}
    extra = {}
    todo = [(p, t) for p in selected for t in p.get("templates", [])]
    while todo:
        p, t = todo.pop()
        if t.get("status") == "merged":
            ref = t["merged_into"]
        elif t.get("status") == "emitted" and t.get("derived"):
            ref = t["derived"]["primary"]
        else:
            continue
        if (ref["record"], ref["key"]) in have or (ref["record"], ref["key"]) in extra:
            continue
        owner = by_record.get(ref["record"])
        mine = [x for x in (owner or {}).get("templates", [])
                if x.get("key") == ref["key"] and x.get("status") == "emitted"]
        if not mine:
            raise EmissionError(f"owner {ref['record'][:16]} of template {t['key'][:16]} "
                                f"is not in the catalog")
        extra[(ref["record"], ref["key"])] = {
            "record_event_id": ref["record"], "family": owner.get("family"),
            "status": "emitted", "carried_for": p.get("record_event_id"),
            "templates": [mine[0]]}
        todo.append((extra[(ref["record"], ref["key"])], mine[0]))
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
                             "refused:locality_witness", "refused:locality_spread",
                             "refused:base_template",
                             "refused:no_base_template", "refused:orientation",
                             "refused:unmapped")


def _molecularity(family, direction) -> int:
    """Reactant count of a family direction (``windows.family_constants``:
    a family whose template is a single tree root, such as R_Recombination,
    lists one template reactant for a bimolecular direction)."""
    from rmgpy.kmc import windows as win
    return dict(win.family_constants(family).directions)[direction]
