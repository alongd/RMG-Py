"""MoM kernel rates from a kMC template catalog (I-213, slice S1 of the
shared-chemistry design).

``kernel_rates(doc, spec, temperature)`` maps every template a
``kmc-template-catalog/1`` loads at ``temperature`` either to a method-of-
moments kernel or to a named kMC-only entry, and gives each kernel's rate
coefficients in the moment unit.  Nothing is fitted and nothing is dropped:
a template the rules below cannot place is listed in ``kmc_only`` with its
reason.

**Kernel type** (from the template's edit graph alone; ``H`` is the only
element read, as the hydrogen-transfer marker; ``delta`` = radical
electrons after minus before):

=========================  ===========  =====  ==========================================
type                       molecularity delta  heavy-atom bond edits; H transfer
=========================  ===========  =====  ==========================================
homolysis                  1            +2     one single bond broken
recombination              2            -2     one single bond formed across participants
disproportionation         2            -2     one pi bond formed; H moves
reverse_disproportionation 2            +2     one pi bond broken; H moves
h_transfer                 2            0      none; H moves
intra_h_transfer           1            0      none; H moves
beta_scission              1            0      one single bond broken, one pi bond formed
addition                   2            0      one single bond formed across, one pi broken
=========================  ===========  =====  ==========================================

**Site classes** (where a participant's edited atoms land on declared
chains).  States are built from the declaration alone (``materialize`` of one
chain of n declared units), by edits on its backbone (``origin.is_backbone``):

* ``chain``: a pristine chain; its count is the per-chain polynomial
  ``s(n) = a*n + b`` in the moment unit, read from the melt engine's own index
  (``MeltEngine._scan``: allowed unimolecular transitions, bimolecular sites per
  branching kind) on chains of every primitive period count up to ``n_max``
  (I-266: a chain of n periods that is no whole number of declared units is
  the first n periods of a declared chain, cut and capped as the chain probe
  below, so the odd-period alkanes of a two-period unit are states too) and
  required to be exactly affine on the trailing lengths (refused otherwise:
  never fitted);
* ``radical_end``: the radical atom left by cutting one backbone single bond
  (every backbone bond, both fragments);
* ``radical_interior``: one H removed from an interior backbone atom;
* ``alkene_terminal``: a cut fragment whose radical atom forms a pi bond
  with its backbone neighbour (that neighbour loses an H);
* ``alkene_internal``: one H removed from each of two bonded interior
  backbone atoms, which then share a pi bond.

Every feature is built at every position of every chain length.  A feature
class is named by the feature atoms' element and hydrogen count
(``radical_interior:CH``; a short fragment gives its own class, e.g. the
methyl radical ``radical_end:CH3``).  Its count is per feature: the sites
whose edited atoms include a feature atom.  ``per_feature`` is the far-field
count (a feature in the middle of the chain: a cut of a middle backbone bond,
or the middle backbone atoms and bonds), required to be the same on the
three longest chains (refused otherwise); ``profile`` lists every position
whose counts differ from it, keyed ``fragment_backbone_atoms=m`` (an end
feature on a fragment of m backbone atoms) or ``chain_backbone_atoms=nb,
end_distance=d`` (an interior feature d backbone bonds from the nearer end
of a chain of nb backbone atoms).  A template that matches only near an end
(a context template of the catalog) has ``per_feature`` 0 and a profile: it
is placed in its kernel with that profile, not dropped.

**Kernel id**: ``<type>:<role>=<class>+...`` with the participants sorted by
role (``roles``: ``donor`` loses the H or forms the pi bond, ``acceptor``
gains it, ``alkene`` opens its pi bond, ``radical`` couples, ``substrate`` for
a unimolecular kernel), so a template written in either participant order
gives one kernel and ``S0``/``S1`` follow the sorted order.

**Moment unit**: the declared backbone's primitive period
(``declared_chain.ChainGeometry``: ``d_rep / period`` periods per declared
unit).  ``n`` in every polynomial counts primitive periods.

**Rates** (SI; ``w`` = the loaded template weight at ``temperature``, per
matched transition, degeneracy removed at emission; ``divisor`` = the melt
engine's ``_Class.divisor``, the edit symmetry times 2 for a role swap):

* unimolecular, ``chain`` class: ``R = sum_T w_T*(a_T*mu1 + b_T*mu0)``,
  reported as ``k_mu1``, ``k_mu0`` (1/s); feature class:
  ``R = sum_T w_T*c_T * mu0(featured pool)``, ``k_per_feature`` (1/s);
* bimolecular: ``R = sum_T (w_T/divisor_T) * S0 * S1 / V`` (mol/s), with
  ``S_p`` the participant's sites in mol (``a*mu1 + b*mu0`` or ``c*mu0`` of
  the pool carrying its class) and V the melt volume (m3); the
  same-component pair term of the melt engine (O(1/N)) is not in this form.

``phi_thinned`` marks a radical-creating template, which the melt engine
accepts with the cage-escape factor; the coefficients here are intrinsic.

``kmc_only`` reasons: ``no_kernel`` (no type above), ``branching`` (every
site pairing, or every transition, leaks: the product is not a linear chain),
``small_species_participant`` (a participant with no port, i.e. a whole
molecule, matches no state built here: a small species that no piece of a
declared chain is), ``no_site_on_declared_states`` (any other participant
with no site on any state built here, e.g. a homolysis at a branch point),
``deferred`` and ``deferred_not_a_product`` (below).

**Deduplicated overlaps** (I-271).  A template that defers to others
(``defers_to``, ``window_matcher.deferred``) makes only the transitions none
of them makes.  A unimolecular class's count is the engine's own (its index
holds only those transitions).  A bimolecular class's sites carry the engine's
marker (``melt_engine._Shared``: the classes it defers to that have the same
site on that state), and a pair of sites whose markers share a class is that
class's transition (``melt_engine._suppressed``).  So a deferring template's
term at a class pairing counts the pairs of marker sets sharing no class:
when those pairs are one product (a set of participant 0's markers with a set
of participant 1's), ``S0``/``S1`` and the laws (``law``, ``law_near_end``,
``law_end_distance``; I-275) are read on those sites only,
and ``sites_deferred_to`` names the classes the term leaves out; when every
pair is suppressed the pairing is no term; a template with no other pairing
is ``deferred``.  When they are no single product the template is
``deferred_not_a_product``: a kernel term is one product of site counts.
The participant counts of a template's row are all its sites.

**Product laws** (``kmc-kernel-rates/2``, I-214; what the kernel references
of ``kernel_reference`` integrate).  Each kernel term carries ``law``: per
participant (``S0``/``S1`` in the sorted order) what every counted site
makes, read by firing the melt engine's own transitions on declared probe
states (``_Probes``: the first n periods of a materialized chain, cut and
capped, with the feature built as above) and classifying each product as a
declared state (pool family, label, n; isomorphic to the probe of that
state, refused otherwise).  A bimolecular participant's product is read with
a partner of the other class and must not depend on the partner's site; a
product holding both whole reactants is a join (``["join", family,
gamma]``, length n0 + n1 + gamma backbone atoms).  ``anchor``: ``chain``
(sites along a pristine chain), ``end`` (an end feature) or ``interior``.
Chain and end laws give ``exact`` outcomes ``[[multiplicity, [[label, n],
...]], ...]`` on every n up to ``exact_to_n`` and a far rule beyond
(``predict``): products as ``[family, alpha, beta, gamma]``, length ``alpha*L
+ beta*k + gamma`` backbone atoms for a reactant of L with the site anchored
at backbone atom k (``bulk``, and the anchors within ``edge`` of either end,
``first``/``last``), required to reproduce the exact outcomes on the
trailing lengths (an interior product compared by family and length, the
resolution of the codes: one length of an interior family can carry several
classes, e.g. a period of a CH2 and a substituted CH).  Interior laws give
the far-field outcomes relative to the feature position.  One length rule
for every participant: a law reads only the probe states carrying its class
(an interior probe's class is read at its position), and needs the class at
the law's reach (``laws_exact_to_n``) and, along a chain or at an end, on
``EXACT_POINTS`` lengths; an end class carried at short lengths only (a
length label of its family, e.g. the methyl radical at 1 period, ethylene at
2) has no longer state, so its law is exact on those lengths
(``fixed_lengths``) with no far outcome.  ``law_inert``: an interior
participant with no far-field site (the rate law above gives 0);
``law_refused``: a product that is no declared state (e.g. a vinyl radical),
a class no probe state carries, or a far field the law cannot state (an
interior law whose neighbouring far-field states differ, a chain law whose
end effects reach the middle of the longest probe, a far rule that does not
reproduce the exact outcomes or gives an impossible product): every such
condition refuses its law by name, never the table.
``families``: per pool family its
label per length, the molar-mass offset (M = n * unit mass + offset), the
element-composition offset (``kmc-kernel-rates/3``, I-217: atoms = n *
``moment_unit.composition`` + ``composition_offset``, exact rationals as
strings, keyed like the melt observables' compositions) and the lengths that
are declared species.  ``overrides``: each depository
override of the term's template at its exact sizes (reactant states with
interior positions, products, matched multiplicity, weight, the catalog's
depository ``entry``); an override
placed in no term is listed in ``overrides_unplaced``.  A site count is the
template's whatever weight a matched size carries.

**Positions** (I-272; optional keys of ``kmc-kernel-rates/3``, every other
key unchanged).  A law any of whose products carries an interior feature
(``radical_interior``, ``alkene_internal``) also has ``positioned``: the same
outcomes, read from the same fired transitions, resolved by where each
feature sits, its product code followed by one of ``{"end_distance": d}``
(d backbone bonds from the nearer end of the product: an edge or end
anchor, or d from the cut end of a fragment), ``{"bulk": [p, q]}`` (a bulk
row: over the bulk anchors the feature's first atom is uniform on the
product's backbone atoms [p, Lp - w - q], Lp atoms, the feature spanning w),
``{"offset": c}`` (an interior law: c backbone atoms past the reactant's own
feature) or ``{"junction": true}`` (a join; where it sits depends on both
reactants).  ``positioned`` has ``exact`` (every probed length) and ``far``
(the shape of the law's own far rule; a chain law its own ``edge``), and its
far rule must reproduce the exact positions on the trailing lengths
(``predict_placed``); otherwise the law says ``positions_refused`` and why.
**End-distance laws** (I-272): a term marked ``law_inert`` keeps ``law``
None and its reason, and gains ``law_end_distance`` (or
``law_end_distance_refused``): per participant (``S0``/``S1``, or the
law itself for one participant) the far-field law of a participant with
one, and for an interior participant with sites only near an end (anchor
``end_distance``) its ``at``: per end distance d with sites on chains of
``from_n`` periods and more, the sites of a feature at d (the count read on
the longest probe with the feature at d, checked against the participant's
profile there; on a chain type that is not its own reverse the key is
``head_distance=d`` or ``tail_distance=d``, below), its positioned outcomes
and, for one participant,
``k_per_feature`` (1/s per feature at d); and ``exact``: on every shorter
length, the same per canonical position (an entry the probes cannot read
says ``refused``).  Nothing is fitted: the sites and products are the
engine's own on declared probe states, the weights the template's.
**Near-end laws** (I-273): a law term whose interior participant has
far-field sites gains ``law_near_end`` (per slot, or the law itself for one
participant; ``law_near_end_refused`` with the reason when the probes
cannot read it) where that participant's sites near an end are not its far
law's: an end-distance law (as above) whose ``at`` lists every end distance
d where the sites of a feature at d, or their positioned outcomes, differ
from the far law at d (zero sites included; scanned one past the reach of
the participant's profile on the longest probe, and refused if they still
differ there), and whose ``exact`` gives every canonical position of the
lengths shorter than where the profile settles.  Elsewhere the far law
holds (``predict_near``).

**Chain types** (I-279, ``kmc-kernel-rates/4``).  A pristine piece is n
primitive periods plus a residue r < nbb backbone atoms, cut from a declared
chain at period offset s and capped by an H at each cut end; (s, r) is its
chain type.  Each piece is told apart by its end radicals alone (no polymer
literal): the piece (s, r) is another piece read backwards when its
start-end radical is isomorphic to that piece's last-end radical, and it is
its own reverse when its two end radicals are isomorphic.  A type is named
by the declaration keys of its two end atoms (``c1~c2``; PP: the declared
``c1~c2``, ``c1~c1`` and ``c2~c2``, the last two their own reverse).  A
declaration with one chain type that is its own reverse (PE) keeps the
``kmc-kernel-rates/3`` state model: none of the keys below appears, and the
table is byte-identical to its ``/3`` table but for ``schema``.  Otherwise:
``chain_types`` gives per type its ``start_offset``, ``residue``,
``symmetric`` (its own reverse), ``declared`` and per family the classes it
carries, offsets and declared species (``_chain_types_doc``); a class label
of a type other than the declared one is ``<class>/<type>`` (the chain pool
``chain/<type>``), with lengths n in periods past the residue; both ends of
every cut are probed, so an end law on a type that is not its own reverse
says ``side: "start"`` when its class is the start end's (the tail-side end
of a cut); interior positions on such a type are oriented, never folded:
a product's position is ``{"head_distance": d}`` or ``{"tail_distance":
d}`` (its feature d backbone bonds from the head or the tail, the nearer
one), and an end-distance law's ``at`` keys are ``head_distance=d`` /
``tail_distance=d``; a type that is its own reverse keeps ``end_distance``,
and its far-field interior states are read off the piece's mirror: a state
where a product holding the whole reactant has its feature past the middle
of the piece (read from the other end than the reactant's feature, e.g. the
two alkenes of a CH donor at the centre atom of a c1~c1 piece, mirror images
both read C=CH) is passed over for the next middlemost (I-279b).
Product codes carry their class (``[label, alpha, beta, gamma]``), so a far
rule is checked by class; a chain law whose bulk depends on the anchor's
place in the period gives ``bulk_phases`` (one row per phase k % nbb, the
anchor k read from the reactant's head) instead of ``bulk``.  A product
whose period is not continued across a junction (head-to-head or
tail-to-tail) is a defect, refused by name.  ``unreached`` lists the
participant classes and the chain types no term's law reads.  A reader of
the ``/3`` state model refuses such a table by name
(``state_model_beyond_v3``).

With a ``kmc-melt-run/1`` manifest, ``hazard`` gives the fired counts per
template, kernel and kMC-only entry (gross), the counts netted over each
reciprocal record pair (``reverse_of``), and the run's pooled leak fraction.

Progress goes to stderr: one line per chain length counted and one per law
term (its template, classes, outcome and seconds).  The table reads the
declaration's chemistry only, never its ensemble.
"""

import argparse
import hashlib
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import window_matcher as wm
from rmgpy.kmc.declared_chain import ChainGeometry
from rmgpy.kmc.melt_state import SetBond, SetRadical, materialize

SCHEMA = "kmc-kernel-rates/4"
EXACT_POINTS = 3          # trailing lengths on which a polynomial / feature count must hold exactly
EXTRA_UNITS = 3           # chain length beyond twice the largest window, in primitive periods


class KernelRatesError(ValueError):
    """The catalog, the declaration or a count is refused."""


# --------------------------------------------------------------------------
# kernel type from the edit graph
# --------------------------------------------------------------------------
def kernel_type(template, cls) -> Optional[str]:
    """The kernel type of one loaded ``window_matcher.Template`` (module
    table), or None."""
    atoms = template.atoms
    heavy = lambda i: atoms[i].element != "H"
    single_formed = single_broken = pi_formed = pi_broken = 0
    h_moves = False
    for (i, j), new in template.bond_final.items():
        old = template.bond.get((i, j), 0.0)
        if not (heavy(i) and heavy(j)):
            h_moves = True
            continue
        if old == 0.0 and new == 1.0:
            single_formed += 1
        elif old == 1.0 and new == 0.0:
            single_broken += 1
        elif old == 1.0 and new == 2.0:
            pi_formed += 1
        elif old == 2.0 and new == 1.0:
            pi_broken += 1
        else:
            return None
    uni, delta = not cls.bimolecular, cls.radical_delta
    edits = (single_formed, single_broken, pi_formed, pi_broken, h_moves)
    table = {
        (True, 2, (0, 1, 0, 0, False)): "homolysis",
        (False, -2, (1, 0, 0, 0, False)): "recombination",
        (False, -2, (0, 0, 1, 0, True)): "disproportionation",
        (False, 2, (0, 0, 0, 1, True)): "reverse_disproportionation",
        (False, 0, (0, 0, 0, 0, True)): "h_transfer",
        (True, 0, (0, 0, 0, 0, True)): "intra_h_transfer",
        (True, 0, (0, 1, 1, 0, False)): "beta_scission",
        (False, 0, (1, 0, 0, 1, False)): "addition",
    }
    found = table.get((uni, delta, edits))
    if found == "recombination" and not cls.cross_bonds:
        return None
    return found


def roles(template, ktype) -> List[str]:
    """Per participant, its role in the kernel: ``donor`` (loses an H, or
    forms the pi bond in a disproportionation), ``acceptor`` (gains the H),
    ``alkene`` (its pi bond opens), ``radical`` (a radical that bonds), or
    ``substrate`` (a unimolecular kernel)."""
    atoms = template.atoms
    parts = template.participants
    if len(parts) == 1:
        return ["substrate"]
    out = {}
    for (i, j), new in template.bond_final.items():
        old = template.bond.get((i, j), 0.0)
        for a, b in ((i, j), (j, i)):
            if atoms[b].element == "H" and atoms[a].element != "H":
                if old and not new:
                    out.setdefault(atoms[a].participant, "donor")
                elif new and not old:
                    out.setdefault(atoms[a].participant, "acceptor")
        if atoms[i].element != "H" and atoms[j].element != "H":
            if old == 1.0 and new == 2.0:
                out[atoms[i].participant] = "donor"
            elif old == 2.0 and new == 1.0:
                out[atoms[i].participant] = "alkene"
    if ktype in ("disproportionation",):
        for p in parts:
            if out.get(p) != "donor":
                out[p] = "acceptor"
    elif ktype in ("recombination", "addition"):
        for p in parts:
            out.setdefault(p, "radical")
    return [out.get(p, "substrate") for p in parts]


# --------------------------------------------------------------------------
# declared states
# --------------------------------------------------------------------------
def _single_chain_spec(spec, units: int):
    from rmgpy.polymer_spec import LinearHomopolymerSpec
    mapping = json.loads(json.dumps(spec.as_mapping()))
    mapping["initial_ensemble"] = {"dp_counts": {units: 1}}
    mapping["initial_mobile_counts"] = {}
    return LinearHomopolymerSpec.from_mapping(mapping)


def _backbone_order(state) -> List[int]:
    """The chain's backbone atoms, from one end to the other."""
    origin = state.origin
    bb = [a for a in state.atoms if origin[a].kind == "chain" and origin[a].is_backbone]
    nbrs = {a: [b for b in state.atoms[a].edges if b in set(bb)] for a in bb}
    ends = sorted(a for a in bb if len(nbrs[a]) <= 1)
    if not bb or (len(bb) > 1 and len(ends) != 2):
        raise KernelRatesError("the declared chain's backbone is not a simple path")
    path, prev = [ends[0]], None
    while len(path) < len(bb):
        nxt = [b for b in nbrs[path[-1]] if b != prev]
        prev = path[-1]
        path.append(nxt[0])
    return path


def _hydrogens(state, a) -> List[int]:
    return sorted(b for b in state.atoms[a].edges if state.atoms[b].element == "H")


def _label(state, atoms) -> str:
    parts = []
    for a in atoms:
        h = sum(1 for b in state.atoms[a].edges if state.atoms[b].element == "H")
        parts.append(state.atoms[a].element + ("H" if h == 1 else f"H{h}" if h else ""))
    return "=".join(parts)


def _remove_h(state, a) -> list:
    h = _hydrogens(state, a)
    if not h:
        return []
    node = state.atoms[a]
    return [SetBond(a, h[0], 1.0, 0.0), SetRadical(a, node.radical, node.radical + 1),
            SetRadical(h[0], 0, 1)]


# --------------------------------------------------------------------------
# counting through the melt engine's own index
# --------------------------------------------------------------------------
def _uni_all(entry) -> Dict[int, Tuple[tuple, int]]:
    """class index -> (allowed transitions, leak count) of an index entry, the
    transitions a depository override re-weights included: a site count is the
    template's, whatever weight a matched size carries."""
    out = {ci: (tuple(allowed), leak) for ci, (allowed, leak) in entry.uni.items()}
    for (ci, _o), (allowed, leak) in sorted(entry.uni_o.items()):
        mine, n = out.get(ci, ((), 0))
        out[ci] = (tuple(sorted(mine + tuple(allowed), key=lambda t: t.edit_key)), n + leak)
    return out


class _Counter:
    def __init__(self, engine):
        self.engine = engine
        self.leak = Counter()        # class index -> unimolecular leak transitions seen anywhere

    def counts(self, state, cid, feature: frozenset = frozenset()) -> Counter:
        """(class index, participant, kind) -> number of sites of component
        ``cid``; with ``feature``, only sites/transitions whose edited atoms
        include a feature atom.  Unimolecular: participant 0, kind None."""
        entry = self.engine._scan([cid], state=state)[cid]
        out = Counter()
        for ci, (allowed, leak) in _uni_all(entry).items():
            self.leak[ci] += leak
            for tr in allowed:
                touched = set()
                for e in tr.edits:
                    touched.update((e.a, e.b) if isinstance(e, SetBond) else (e.atom,))
                if not feature or touched & feature:
                    out[(ci, 0, None)] += 1
        for (ci, p, kind), sites in entry.sites.items():
            n = sum(1 for s in sites if not feature or set(s) & feature)
            if n:
                # a depository override's tag is a weight, not a site class; a deferral marker is kept
                out[(ci, p, kind[:2] + (kind[-1],) if isinstance(kind[-1], me._Shared) else kind[:2])] += n
        return out


def _marker(kind) -> frozenset:
    """The classes a counted site kind defers its shared pairs to (empty for a
    site of no deferring class, or a unimolecular count)."""
    return frozenset(kind[-1]) if kind and isinstance(kind[-1], me._Shared) else frozenset()


def _unshared_blocks(m0, m1) -> List[Tuple[frozenset, frozenset]]:
    """The pairs of site markers (participant 0's, participant 1's) a class
    fires itself, those sharing no class (``melt_engine._suppressed``), as
    products: [(markers of participant 0, markers of participant 1)]."""
    by = defaultdict(set)
    for a in m0:
        allowed = frozenset(b for b in m1 if not a & b)
        if allowed:
            by[allowed].add(a)
    return sorted(((frozenset(a), b) for b, a in by.items()), key=lambda x: sorted(map(sorted, x[0])))


def _affine(series: Dict[int, int]) -> dict:
    """a, b of s(n) = a*n + b exact on the trailing EXACT_POINTS lengths, and
    the smallest length from which it holds; refused when not exact."""
    pts = sorted(series.items())
    if len(pts) < EXACT_POINTS:
        raise KernelRatesError("too few chain lengths to establish a site polynomial")
    (n1, y1), (n2, y2) = pts[-2], pts[-1]
    a = (y2 - y1) / (n2 - n1)
    b = y1 - a * n1
    exact = [n for n, y in pts if abs(y - (a * n + b)) <= 1e-9]
    tail = [n for n, _y in pts[-EXACT_POINTS:]]
    if any(n not in exact for n in tail):
        raise KernelRatesError(f"site counts {dict(pts)} are not affine in n on the trailing lengths")
    start = tail[0]
    for n, _y in reversed(pts[:-EXACT_POINTS]):
        if n in exact:
            start = n
        else:
            break
    return {"a": a, "b": b, "exact_from_n": start}


def _chain_state(spec, n: int, ppu: int, nbb: int):
    """(state, component) of a pristine chain of n primitive periods: the
    declared chain of n / ppu units when n is a whole number of them, else
    the first n periods of a longer declared chain cut from the rest, the cut
    end capped by an H of the rest (``_Probes.put``'s chain probe)."""
    if n % ppu == 0:
        state = materialize(_single_chain_spec(spec, n // ppu))
        return state, state.components.component_of(_backbone_order(state)[0])
    state = materialize(_single_chain_spec(spec, n // ppu + 2))
    path = _backbone_order(state)
    L = n * nbb
    u, v, q = path[L - 1], path[L], path[L + 1]
    hq = _hydrogens(state, q)
    if state.atoms[u].edges.get(v) != 1.0 or state.atoms[v].edges.get(q) != 1.0 or not hq:
        raise KernelRatesError(f"no pristine chain of {n} periods: the declared chain is not cut and capped there")
    state.apply([SetBond(u, v, 1.0, 0.0), SetBond(q, hq[0], 1.0, 0.0), SetBond(u, hq[0], 0.0, 1.0),
                 SetBond(v, q, 1.0, 2.0)])
    return state, state.components.component_of(path[0])


def _piece(spec, s: int, r: int, n: int, ppu: int, nbb: int):
    """(state, backbone) of a pristine piece of chain type (s, r) (module doc,
    chain types) of n periods: ``n * nbb + r`` backbone atoms from period
    offset s of a longer declared chain, each cut end capped by an H of the
    rest (``_Probes.put``'s chain probe); a piece at offset 0 keeps the
    declared head."""
    L = n * nbb + r
    i0 = _start_index(s, nbb) if s else 0
    state = materialize(_single_chain_spec(spec, -(-(i0 + L + 2) // (ppu * nbb)) + 1))
    path = _backbone_order(state)
    for a, b, c in ([(path[i0], path[i0 - 1], path[i0 - 2])] if i0 else []) + \
            [(path[i0 + L - 1], path[i0 + L], path[i0 + L + 1])]:
        hc = _hydrogens(state, c)
        if state.atoms[a].edges.get(b) != 1.0 or state.atoms[b].edges.get(c) != 1.0 or not hc:
            raise KernelRatesError(f"no pristine piece of {L} backbone atoms at period offset {s}: the declared "
                                   "chain is not cut and capped there")
        state.apply([SetBond(a, b, 1.0, 0.0), SetBond(c, hc[0], 1.0, 0.0), SetBond(a, hc[0], 0.0, 1.0),
                     SetBond(b, c, 1.0, 2.0)])
    return state, path[i0:i0 + L]


def _states(spec, n_values, counter, ppu: int, nbb: int, probes=None):
    """``chain``: {class: {(class index, participant, kind): {n: count}}} on
    pristine pieces of every chain type (class ``chain``, or ``chain/<type>``
    for a type other than the declared one; ``probes``' types) of every n
    primitive periods up to the longest declared chain (``_chain_state``,
    ``_piece``); ``features``: {class label: {"far": {units: counts},
    "profile": {position: [counts, ...]}}} (``_features``, on the pieces of
    every chain type of ``n_values`` units)."""
    chains: Dict[str, Dict[Tuple, Dict[int, int]]] = {}
    features: Dict[str, dict] = {}
    longest = max(n_values)
    types = [None] + ([] if probes is None else sorted(t for t in probes.sr if t != probes.declared))
    for tau in types:
        chain = chains.setdefault("chain" if tau is None else f"chain/{tau}", defaultdict(dict))
        s, r = (0, 0) if tau is None else probes.sr[tau]
        for n in range(0 if r else 1, longest * ppu + 1):
            _progress(f"site counts on {n}/{longest * ppu} periods" + ("" if tau is None else f" of type {tau}"))
            if tau is None:
                state, cid = _chain_state(spec, n, ppu, nbb)
            else:
                state, seg = _piece(spec, s, r, n, ppu, nbb)
                cid = state.components.component_of(seg[0])
            for key, v in counter.counts(state, cid).items():
                chain[key][n] = v
            for key in chain:
                chain[key].setdefault(n, 0)
    for tau in types:
        for units in n_values:
            _progress(f"feature counts on {units}/{longest} declared units" + ("" if tau is None else f" of type {tau}"))
            _features(spec, units, counter, features, far_ok=(units > longest - EXACT_POINTS), probes=probes, tau=tau)
    return chains, features


def _features(spec, units, counter, features, far_ok, probes=None, tau=None):
    """Feature counts on one chain of ``units`` declared units (a piece of
    chain type ``tau`` of as many periods), each feature on a fresh state, at
    every position: ``("fragment_backbone_atoms", m)`` for an end feature on a
    fragment of m backbone atoms, and ``("chain_backbone_atoms", nb,
    "end_distance", d)`` for an interior feature d backbone bonds from the
    nearer end of a chain of nb backbone atoms (``head_distance`` /
    ``tail_distance`` on a chain type that is not its own reverse).  An end
    feature's class is that of its fragment's chain type.  With ``far_ok`` the
    features in the middle of the chain (a cut of a middle bond; the middle
    backbone atoms or bonds) also give the far-field counts."""
    legacy = probes is None or probes.legacy
    if tau is None:
        def fresh():
            s = materialize(_single_chain_spec(spec, units))
            return s, _backbone_order(s)
        s0 = 0
    else:
        s0, r0 = probes.sr[tau]

        def fresh():
            return _piece(spec, s0, r0, units * probes.ppu, probes.ppu, probes.nbb)
    nb = len(fresh()[1])
    per_unit = max(1, nb // units)
    middle = range((units // 2) * per_unit, (units // 2) * per_unit + per_unit)
    sym = legacy or probes.symmetric(tau)
    if sym and not legacy:
        # a type that is its own reverse reads a feature from its nearer end: the middle unit is the last whole
        # unit on the head side of the centre, so every class of the unit is read the same way at every length
        start = max(0, (nb - 2) // 2 - per_unit + 1)
        middle = range(start, start + per_unit)

    def labelled(kind, atoms, s, r):
        """The class label of a feature on a piece of period offset s and residue r."""
        if legacy:
            return f"{kind}:{atoms}"
        name = probes._tau(probes.canon[(s % probes.nbb, r % probes.nbb)][0])
        return f"{kind}:{atoms}" + ("" if name is None else f"/{name}")

    def note(name, counts, position, far):
        slot = features.setdefault(name, {"far": {}, "profile": {}})
        if far and far_ok:
            slot["far"].setdefault(units, counts)
        slot["profile"].setdefault(position, []).append(counts)

    cuts = range(1, nb)
    for k in cuts:
        for side in (0, 1):
            s, p = fresh()
            u, v = p[k - 1], p[k]
            if s.atoms[u].edges.get(v) != 1.0:
                continue
            s.apply([SetBond(u, v, 1.0, 0.0), SetRadical(u, s.atoms[u].radical, s.atoms[u].radical + 1),
                     SetRadical(v, s.atoms[v].radical, s.atoms[v].radical + 1)])
            if side == 0:
                r, q, m, start = u, (p[k - 2] if k >= 2 else None), k, s0
            else:
                r, q, m, start = v, (p[k + 1] if k + 1 < nb else None), nb - k, s0 + k
            far = (k - 1) in middle
            position = ("fragment_backbone_atoms", m)
            note(labelled("radical_end", _label(s, [r]), start, m),
                 counter.counts(s, s.components.component_of(r), frozenset([r])), position, far)
            if q is not None and _hydrogens(s, q):
                hq = _hydrogens(s, q)[0]
                s.apply([SetBond(q, hq, 1.0, 0.0), SetRadical(hq, 0, 1), SetBond(r, q, 1.0, 2.0),
                         SetRadical(r, s.atoms[r].radical, s.atoms[r].radical - 1)])
                note(labelled("alkene_terminal", _label(s, [r, q]), start, m),
                     counter.counts(s, s.components.component_of(r), frozenset([r, q])), position, far)

    for j in range(1, nb - 1):
        far = j in middle
        if sym:
            position = alkene_position = ("chain_backbone_atoms", nb, "end_distance", min(j, nb - 1 - j))
        else:
            position = ("chain_backbone_atoms", nb) + _oriented_pos(j, nb, 1)
            alkene_position = ("chain_backbone_atoms", nb) + _oriented_pos(j, nb, 2)
        s, p = fresh()
        edits = _remove_h(s, p[j])
        if edits:
            s.apply(edits)
            note(labelled("radical_interior", _label(s, [p[j]]), s0, nb),
                 counter.counts(s, s.components.component_of(p[j]), frozenset([p[j]])), position, far)
        if j + 1 < nb - 1:
            s, p = fresh()
            a, b = p[j], p[j + 1]
            ha, hb = _hydrogens(s, a), _hydrogens(s, b)
            if ha and hb and s.atoms[a].edges.get(b) == 1.0:
                s.apply([SetBond(a, ha[0], 1.0, 0.0), SetBond(b, hb[0], 1.0, 0.0),
                         SetRadical(ha[0], 0, 1), SetRadical(hb[0], 0, 1), SetBond(a, b, 1.0, 2.0)])
                pair = [b, a] if not legacy and sym and j > nb - 2 - j else [a, b]
                note(labelled("alkene_internal", _label(s, pair), s0, nb),
                     counter.counts(s, s.components.component_of(a), frozenset([a, b])), alkene_position, far)


# --------------------------------------------------------------------------
# product laws (I-214): what every counted site makes
# --------------------------------------------------------------------------
END_FAMILIES = ("radical_end", "alkene_terminal")
INTERIOR_FAMILIES = ("radical_interior", "alkene_internal")
FAMILIES = ("chain",) + END_FAMILIES + INTERIOR_FAMILIES


class _OffDeclared(KernelRatesError):
    """A product (or a probe) that is no declared state."""


def _declared_counts(composition) -> Counter:
    """Atom counts of a declared composition by (element, isotope or None)."""
    return Counter({(symbol, None if isotope == -1 else isotope): count
                    for symbol, isotope, count in composition.entries})


def family_of(label: str) -> str:
    return label.split("/")[0].split(":")[0]


def type_of(label: str) -> Optional[str]:
    """The chain type a class label names (``<class>/<type>``), None for the declared type."""
    return label.split("/", 1)[1] if "/" in label else None


def _type_entry(token: str, types: dict) -> dict:
    """The ``chain_types`` entry of a class label (or product code token)."""
    name = type_of(token)
    if name is None:
        name = next(k for k, v in types.items() if v["declared"])
    return types[name]


def _residue(token: str, types: Optional[dict]) -> int:
    return 0 if types is None else _type_entry(token, types)["residue"]


def _symmetric(token: str, types: Optional[dict]) -> bool:
    return True if types is None else _type_entry(token, types)["symmetric"]


def _start_index(s: int, nbb: int) -> int:
    """The first backbone index of a declared chain at period offset s with two backbone atoms before it (a cut
    there is capped by an H of the atom two back)."""
    return s + nbb * max(0, -(-(2 - s) // nbb))


def _oriented_pos(x: int, lp: int, w: int) -> Tuple[str, int]:
    """An oriented position (the feature's first atom x atoms from the head of a product of lp backbone atoms, the
    feature spanning w) by its nearer end: ``("head_distance", x)`` or ``("tail_distance", lp - w - x)``."""
    return ("head_distance", x) if x <= lp - w - x else ("tail_distance", lp - w - x)


def _product_pos(x: int, lp: int, w: int, symmetric: bool) -> Tuple[str, int]:
    """A product's position code: folded to the nearer end on a chain type that is its own reverse, else oriented."""
    return ("end_distance", min(x, lp - w - x)) if symmetric else _oriented_pos(x, lp, w)


def state_model_beyond_v3(table: dict) -> Optional[str]:
    """Why a kernel-rates table cannot be read in the ``kmc-kernel-rates/3`` state model (one chain type, its own
    reverse, codes by family), or None when it can: a reader of that model refuses the table by this reason."""
    types = table.get("chain_types")
    if not types:
        return None
    return (f"the table's state model has {len(types)} chain types ({', '.join(sorted(types))}), oriented "
            "positions and product codes by class (kmc-kernel-rates/4); this reader reads one chain type that is "
            "its own reverse (kmc-kernel-rates/3)")


def _chains_spec(spec, units: int, chains: int):
    from rmgpy.polymer_spec import LinearHomopolymerSpec
    mapping = json.loads(json.dumps(spec.as_mapping()))
    mapping["initial_ensemble"] = {"dp_counts": {units: chains}}
    mapping["initial_mobile_counts"] = {}
    return LinearHomopolymerSpec.from_mapping(mapping)


def _walk(state, bb) -> List[int]:
    """The atoms ``bb`` as a simple path, from its smaller end; refused when
    they are not one."""
    bb = set(bb)
    nbrs = {a: [b for b in state.atoms[a].edges if b in bb] for a in bb}
    ends = sorted(a for a in bb if len(nbrs[a]) <= 1)
    if not bb or (len(bb) > 1 and len(ends) != 2) or any(len(v) > 2 for v in nbrs.values()):
        raise _OffDeclared("a backbone that is not a simple path")
    path, prev = [ends[0]], None
    while len(path) < len(bb):
        nxt = [b for b in nbrs[path[-1]] if b != prev]
        if not nxt:
            raise _OffDeclared("a backbone that is not connected")
        prev = path[-1]
        path.append(nxt[0])
    return path


def _paths(state) -> List[List[int]]:
    """The backbone path of every component of a materialized melt, in
    component order."""
    origin = state.origin
    return [_walk(state, [a for a in state.components.atoms(cid)
                          if origin[a].kind == "chain" and origin[a].is_backbone])
            for cid in state.components.ids()]


def _hydrogen_folded(g):
    """``g`` (``melt_engine._graph``) with every hydrogen leaf (an H whose one
    neighbour is not an H) folded into its neighbour's label, as the sorted
    multiset of (bond order, leaf label).  The leaves of one atom are
    interchangeable, so two graphs are isomorphic iff their folds are; VF2 on
    the unfolded graphs branches over them (exponentially in the CH2 and CH3
    groups of a long probe when a form does not match)."""
    import networkx as nx
    leaves = {a for a, d in g.nodes(data=True) if d["label"][0] == "H" and g.degree(a) == 1
              and g.nodes[next(iter(g[a]))]["label"][0] != "H"}
    out = nx.Graph()
    for a, d in g.nodes(data=True):
        if a not in leaves:
            held = sorted(((g[a][b]["order"], g.nodes[b]["label"]) for b in g[a] if b in leaves), key=repr)
            out.add_node(a, label=(d["label"], tuple(held)))
    out.add_edges_from((a, b, {"order": o}) for a, b, o in g.edges(data="order")
                       if a not in leaves and b not in leaves)
    return out


def _isomorphic(g, h) -> bool:
    """``melt_engine._isomorphic`` on the hydrogen-folded graphs."""
    return me._isomorphic(_hydrogen_folded(g), _hydrogen_folded(h))


def _heavy_touched(tr) -> set:
    touched = set()
    for e in tr.edits:
        touched.update((e.a, e.b) if isinstance(e, SetBond) else (e.atom,))
    return touched


class _Probes:
    """Declared probe states of the product laws.  A probe is a piece of a
    materialized chain cut from the rest: ``chain`` (a cut end takes an H of
    the rest, whose cut atoms close a pi bond), ``radical_end`` and
    ``alkene_terminal`` (the feature at a cut end, as ``_features`` builds
    it: at the piece's last backbone atom, or with ``x == "start"`` at its
    first), ``radical_interior`` / ``alkene_internal`` (on the chain piece, at
    backbone atom / bond ``x`` from the piece's head).  A piece of chain type
    ``tau`` (module doc, chain types; None: the declared type) holds
    ``n * nbb + residue`` backbone atoms from its type's start offset: a piece
    of the declared type with its feature at the last atom is the first n
    periods (the declared head end kept); any other piece is cut at both ends.
    Every probe sits on chain 0 of a two-chain melt of ``units`` declared
    units, so the atom ids of one probe are the same in every state that holds
    it (chain 1 carries a partner for bimolecular sites).  A chain probe at a
    whole number of declared units must be the declared chain (isomorphism),
    or no chain probe is declared."""

    REF_PERIODS = 4         # periods of the pieces whose ends tell the chain types apart

    def __init__(self, spec, engine, ppu: int, nbb: int, units: int):
        self.spec, self.engine, self.ppu, self.nbb, self.units = spec, engine, ppu, nbb, units
        self.unit_mass = Fraction(spec.repeat_unit.composition.molar_mass_g_mol()) / ppu
        self._forms: Dict[Tuple, Tuple] = {}
        self._verified = set()
        self.labels: Dict[Tuple[str, int], str] = {}
        self.labels_t: Dict[Tuple, str] = {}
        self.offsets: Dict[str, set] = defaultdict(set)
        self.type_offsets: Dict[Tuple[str, str], set] = defaultdict(set)
        self.composition = {f"{el}{iso or ''}": Fraction(k, ppu)
                            for (el, iso), k in _declared_counts(spec.repeat_unit.composition).items()}
        self.composition_offsets: Dict[str, set] = defaultdict(set)
        self.type_composition_offsets: Dict[Tuple[str, str], set] = defaultdict(set)
        self.species: Dict[Tuple[str, int], str] = {}
        self._entries: Dict[Tuple, Any] = {}
        self._carried: Dict[Tuple, str] = {}
        self.backbone_keys = tuple(ChainGeometry.of(spec).backbone_path)
        self.sr: Dict[str, Tuple[int, int]] = {}
        self.sym: Dict[str, bool] = {}
        self.declared, self.legacy = None, True
        declared = materialize(_single_chain_spec(spec, 1))
        dcid = declared.components.ids()[0]
        key, forms = self.forms("chain", ppu)
        g = me._graph(me._first_rows(declared, dcid))
        if me._comp_key(declared, dcid) != key or not any(_isomorphic(g, h) for h in forms):
            raise KernelRatesError("a cut chain capped by an H of its neighbour is not the declared chain: "
                                   "the declared ends are not hydrogen")
        state, paths = self.melt()
        if len(self.backbone_keys) != nbb or state.origin[paths[0][0]].key != self.backbone_keys[0]:
            raise KernelRatesError("the probe chains' backbone does not start at the declared period's head")
        self.where = {a: (c, i) for c, path in enumerate(paths) for i, a in enumerate(path)}
        self._chain_types()

    # -- chain types --------------------------------------------------------
    def _raw_forms(self, family, s, L, side):
        state, paths = self.melt()
        seg, _feature = self._put(state, paths[0], family, self.REF_PERIODS, s, L, None, side)
        cid = state.components.component_of(seg[0])
        return (me._comp_key(state, cid),
                [me._graph(wm.form_rows(state, cid, f)) for f in wm.component_forms(state, cid)])

    def _chain_types(self):
        """Tell the chain types apart (module doc, chain types): the piece at period offset s of residue r is
        another piece read backwards when its start-end radical is that piece's last-end radical, and its own
        reverse when its two end radicals are one state."""
        nbb, keys = self.nbb, self.backbone_keys
        same = lambda a, b: a[0] == b[0] and any(_isomorphic(a[1][0], h) for h in b[1])
        self.canon: Dict[Tuple[int, int], Tuple[str, bool]] = {}
        for r in range(nbb):
            for s in range(nbb):
                L = self.REF_PERIODS * nbb + r
                start, stop = self._raw_forms("radical_end", s, L, "start"), self._raw_forms("radical_end", s, L, None)
                hit = None
                for name, (s2, r2) in self.sr.items():
                    if r2 != r:
                        continue
                    other = self._raw_forms("radical_end", s2, L, None)
                    if same(start, other):
                        hit = (name, True)
                    elif same(stop, other):
                        hit = (name, False)
                    if hit:
                        break
                if hit is None:
                    name = f"{keys[s]}~{keys[(s + r - 1) % nbb]}"
                    self.sr[name], self.sym[name] = (s, r), same(start, stop)
                    hit = (name, False)
                self.canon[(s, r)] = hit
        self.declared = self.canon[(0, 0)][0]
        self.legacy = len(self.sr) == 1 and self.sym[self.declared]
        for name in self.sr:
            if self.sym[name]:
                continue
            for family in END_FAMILIES:
                try:
                    ends = [self.label(family, self.REF_PERIODS, side, self._tau(name)) for side in (None, "start")]
                except _OffDeclared:
                    continue
                if ends[0] == ends[1]:
                    raise KernelRatesError(f"chain type {name}: both ends carry the class {ends[0]}, so a class "
                                           "label cannot tell them apart")

    def _tau(self, tau) -> Optional[str]:
        """None for the declared chain type, else the type's name."""
        if tau is None or tau == self.declared:
            return None
        if tau not in self.sr:
            raise KernelRatesError(f"no chain type {tau!r} (types: {sorted(self.sr)})")
        return tau

    def sr_of(self, tau) -> Tuple[int, int]:
        tau = self._tau(tau)
        return (0, 0) if tau is None else self.sr[tau]

    def symmetric(self, tau) -> bool:
        tau = self._tau(tau)
        return self.sym[self.declared if tau is None else tau] if self.declared else True

    def types_doc(self) -> Optional[dict]:
        """The chain types as the module functions read them (None in the ``kmc-kernel-rates/3`` state model)."""
        if self.legacy:
            return None
        return {name: {"start_offset": s, "residue": r, "symmetric": self.sym[name], "declared": name == self.declared}
                for name, (s, r) in sorted(self.sr.items())}

    def make_label(self, family, state, feature, tau=None) -> str:
        base = family if family == "chain" else f"{family}:{_label(state, feature)}"
        tau = self._tau(tau)
        return base if tau is None else f"{base}/{tau}"

    # -- building -----------------------------------------------------------
    def melt(self):
        state = materialize(_chains_spec(self.spec, self.units, 2))
        return state, _paths(state)

    def put(self, state, path, family, n, x=None, tau=None) -> Tuple[List[int], Tuple[int, ...]]:
        """Edits making (family, n[, x]) of chain type ``tau`` on ``path``;
        (its backbone, its feature atoms)."""
        tau = self._tau(tau)
        s, r = self.sr_of(tau)
        side = x if family in END_FAMILIES else None
        if side not in (None, "start"):
            raise KernelRatesError(f"an end feature at {side!r}: the ends are the last atom (None) and 'start'")
        return self._put(state, path, family, n, s, n * self.nbb + r, x if family in INTERIOR_FAMILIES else None,
                         side)

    def _put(self, state, path, family, n, s, L, x, side):
        start_cut = s != 0 or side == "start"
        i0 = _start_index(s, self.nbb) if start_cut else 0
        if L < 1 or i0 + L + 2 > len(path):
            raise _OffDeclared(f"no {family} probe of {n} periods on a chain of {len(path)} backbone atoms")
        u, v, q = path[i0 + L - 1], path[i0 + L], path[i0 + L + 1]
        seg = list(path[i0:i0 + L])
        if state.atoms[u].edges.get(v) != 1.0 or state.atoms[v].edges.get(q) != 1.0:
            raise _OffDeclared(f"no single bond to cut after {n} periods")
        rad = lambda a, d: SetRadical(a, state.atoms[a].radical, state.atoms[a].radical + d)
        if start_cut:
            a, b, c = path[i0], path[i0 - 1], path[i0 - 2]
            if state.atoms[a].edges.get(b) != 1.0 or state.atoms[b].edges.get(c) != 1.0:
                raise _OffDeclared(f"no single bond to cut before {n} periods")
            if side != "start":
                hc = _hydrogens(state, c)
                if not hc:
                    raise _OffDeclared("no hydrogen to cap a cut chain")
                state.apply([SetBond(a, b, 1.0, 0.0), SetBond(c, hc[0], 1.0, 0.0), SetBond(a, hc[0], 0.0, 1.0),
                             SetBond(b, c, 1.0, 2.0)])
        if family in END_FAMILIES and side is None:
            edits = [SetBond(u, v, 1.0, 0.0), rad(v, 1)]
            if family == "radical_end":
                state.apply(edits + [rad(u, 1)])
                return seg, (u,)
            w = seg[-2] if L >= 2 else None
            hw = _hydrogens(state, w) if w is not None else []
            if not hw or state.atoms[u].edges.get(w) != 1.0:
                raise _OffDeclared(f"no terminal pi bond on {n} periods")
            state.apply(edits + [SetBond(w, hw[0], 1.0, 0.0), SetRadical(hw[0], 0, 1), SetBond(u, w, 1.0, 2.0)])
            return seg, (u, w)
        hq = _hydrogens(state, q)
        if not hq:
            raise _OffDeclared("no hydrogen to cap a cut chain")
        state.apply([SetBond(u, v, 1.0, 0.0), SetBond(q, hq[0], 1.0, 0.0), SetBond(u, hq[0], 0.0, 1.0),
                     SetBond(v, q, 1.0, 2.0)])
        if family in END_FAMILIES:
            edits = [SetBond(a, b, 1.0, 0.0), rad(b, 1)]
            if family == "radical_end":
                state.apply(edits + [rad(a, 1)])
                return seg, (a,)
            w = seg[1] if L >= 2 else None
            hw = _hydrogens(state, w) if w is not None else []
            if not hw or state.atoms[a].edges.get(w) != 1.0:
                raise _OffDeclared(f"no terminal pi bond on {n} periods")
            state.apply(edits + [SetBond(w, hw[0], 1.0, 0.0), SetRadical(hw[0], 0, 1), SetBond(a, w, 1.0, 2.0)])
            return seg, (a, w)
        if family == "chain":
            return seg, ()
        if family == "radical_interior":
            if x is None or not 1 <= x <= L - 2 or not _hydrogens(state, seg[x]):
                raise _OffDeclared(f"no interior radical at {x} on {n} periods")
            state.apply(_remove_h(state, seg[x]))
            return seg, (seg[x],)
        if family == "alkene_internal":
            if x is None or not (1 <= x and x + 1 <= L - 2):
                raise _OffDeclared(f"no internal pi bond at {x} on {n} periods")
            a, b = seg[x], seg[x + 1]
            ha, hb = _hydrogens(state, a), _hydrogens(state, b)
            if not (ha and hb) or state.atoms[a].edges.get(b) != 1.0:
                raise _OffDeclared(f"no internal pi bond at {x} on {n} periods")
            state.apply([SetBond(a, ha[0], 1.0, 0.0), SetBond(b, hb[0], 1.0, 0.0),
                         SetBond(ha[0], hb[0], 0.0, 1.0), SetBond(a, b, 1.0, 2.0)])
            return seg, (a, b)
        raise KernelRatesError(f"unknown family {family}")

    def positions(self, family, n, tau=None) -> List:
        """The canonical positions of a family on n periods of chain type ``tau``: interior backbone atoms / bonds
        from the head (on a type that is its own reverse only those nearer the head), the ends of an end family
        (None: the last atom; ``"start"``: the first, on a type that is not its own reverse)."""
        tau = self._tau(tau)
        L = n * self.nbb + self.sr_of(tau)[1]
        sym = self.symmetric(tau)
        if family == "radical_interior":
            return [x for x in range(1, L - 1) if not sym or x <= L - 1 - x]
        if family == "alkene_internal":
            return [x for x in range(1, L - 2) if not sym or x <= L - 2 - x]
        if family in END_FAMILIES and not sym:
            return [None, "start"]
        return [None]

    def entry(self, family, n, x=None, chain=0, tau=None):
        """(seg, feature, index entry) of a probe on chain ``chain`` (the
        other chain left whole); the scan is kept."""
        tau = self._tau(tau)
        key = (family, n, x, chain) if tau is None else (family, n, x, chain, tau)
        if key not in self._entries:
            state, paths = self.melt()
            seg, feature = self.put(state, paths[chain], family, n, x, tau)
            cid = state.components.component_of(seg[0])
            self._note(state, cid, family, n, feature, tau, x if family in END_FAMILIES else None)
            self._entries[key] = (seg, feature, self.engine._scan([cid], state=state)[cid])
        seg, feature, entry = self._entries[key]
        return seg, feature, entry

    def _note(self, state, cid, family, n, feature, tau=None, side=None):
        label = self.make_label(family, state, feature, tau)
        if tau is None and side is None:
            if self.legacy or family not in INTERIOR_FAMILIES:
                old = self.labels.setdefault((family, n), label)
                if old != label and family not in INTERIOR_FAMILIES:
                    raise KernelRatesError(f"{family} on {n} periods has two labels ({old}, {label})")
        elif family not in INTERIOR_FAMILIES:
            old = self.labels_t.setdefault((family, tau, n, side), label)
            if old != label:
                raise KernelRatesError(f"{family} of chain type {tau or self.declared} on {n} periods "
                                       f"({side or 'last'} end) has two labels ({old}, {label})")
        offsets = self.offsets[family] if tau is None else self.type_offsets[(family, tau)]
        offsets.add(Fraction(state.components.molar_mass_g_mol(cid)) - self.unit_mass * n)
        atoms = {f"{el}{iso or ''}": k for (el, iso), k in state.components.composition(cid).items()}
        compositions = self.composition_offsets[family] if tau is None else \
            self.type_composition_offsets[(family, tau)]
        compositions.add(tuple(sorted(
            (el, v) for el in set(atoms) | set(self.composition)
            if (v := atoms.get(el, 0) - n * self.composition.get(el, 0)))))
        name = self.engine.identity.identify(state, cid)
        if name is not None:
            self.species[(label, n)] = name

    def forms(self, family, n, x=None, tau=None):
        tau = self._tau(tau)
        key = (family, n, x) if tau is None else (family, n, x, tau)
        if key not in self._forms:
            state, paths = self.melt()
            seg, feature = self.put(state, paths[0], family, n, x, tau)
            cid = state.components.component_of(seg[0])
            self._note(state, cid, family, n, feature, tau, x if family in END_FAMILIES else None)
            self._forms[key] = (me._comp_key(state, cid),
                                [me._graph(wm.form_rows(state, cid, f)) for f in wm.component_forms(state, cid)])
        return self._forms[key]

    def label(self, family, n, x=None, tau=None):
        """The class label of (family, n) of chain type ``tau``: an end family's at its end ``x``, an interior
        family's at its first canonical position."""
        tau = self._tau(tau)
        if family in INTERIOR_FAMILIES and not (tau is None and self.legacy):
            xs = self.positions(family, n, tau)
            if not xs:
                raise _OffDeclared(f"no {family} position on {n} periods")
            return self._label_at(family, n, xs[0], tau)
        side = x if family in END_FAMILIES else None
        if tau is None and side is None:
            if (family, n) not in self.labels:
                xs = self.positions(family, n)
                if not xs:
                    raise _OffDeclared(f"no {family} position on {n} periods")
                self.forms(family, n, xs[0])
            return self.labels[(family, n)]
        key = (family, tau, n, side)
        if key not in self.labels_t:
            self.forms(family, n, x, tau)
        return self.labels_t[key]

    def _label_at(self, family, n, x, tau=None):
        key = (family, n, x) if tau is None else (family, n, x, tau)
        if key not in self._carried:
            state, paths = self.melt()
            seg, feature = self.put(state, paths[0], family, n, x, tau)
            if not self.legacy:
                # the probe's molar mass, composition and declared species (``label`` reads an interior family
                # here only, off the legacy model)
                self._note(state, state.components.component_of(seg[0]), family, n, feature, tau)
            self._carried[key] = self.make_label(family, state, feature, tau)
        return self._carried[key]

    def carries(self, label, n, x=None) -> bool:
        """Whether the probe (family of ``label``, n[, x]) of the label's chain
        type carries the class ``label``: an interior probe's class is read at
        its own position (the positions of one probe can carry different
        classes), an end probe's at its end ``x``."""
        family = family_of(label)
        tau = self._tau(type_of(label))
        try:
            if family not in INTERIOR_FAMILIES:
                return self.label(family, n, x, tau) == label
            return self._label_at(family, n, x, tau) == label
        except _OffDeclared:
            return False

    # -- reading products ---------------------------------------------------
    def _read(self, path) -> List[int]:
        """A product's backbone read head to tail along the probe chains: each chain's atoms in increasing index
        and, across a junction of two chains, the period continued; refused otherwise (a defect)."""
        where = [self.where.get(a) for a in path]
        if None in where:
            raise _OffDeclared("a product backbone atom off the probe chains")
        steps = {j[1] > i[1] for i, j in zip(where, where[1:]) if i[0] == j[0]}
        if len(steps) > 1:
            raise _OffDeclared("a product joined head to head or tail to tail (a defect)")
        if steps == {False}:
            path, where = path[::-1], where[::-1]
        if any(i[0] != j[0] and (j[1] - i[1] - 1) % self.nbb for i, j in zip(where, where[1:])):
            raise _OffDeclared("a product joined head to head or tail to tail (a defect)")
        return list(path)

    def classify(self, state, cid, segs) -> Tuple:
        """(family, label, n, ends, place, backbone atoms, code token) of a
        component made from the probe backbones ``segs``: ``ends`` per
        segment, (holds its first atom, holds its last atom, holds any of
        it); ``place`` where an interior feature sits (``_place``; None for
        another family); the token is the family in the
        ``kmc-kernel-rates/3`` state model, else the class label.  The
        component is read head to tail, its chain type from its first atom's
        period offset and its length's residue; it must be the declared probe
        of its family, chain type, length and position (isomorphism, once per
        backbone signature)."""
        atoms = set(state.components.atoms(cid))
        path = _walk(state, [a for seg in segs for a in seg if a in atoms])
        L = len(path)
        flipped = False
        if self.legacy:
            if L % self.nbb:
                raise _OffDeclared(f"a product of {L} backbone atoms is off the period grid")
            tau = None
        else:
            path = self._read(path)
            name, flipped = self.canon[(self.where[path[0]][1] % self.nbb, L % self.nbb)]
            if flipped:
                path = path[::-1]
            tau = self._tau(name)
        sym = self.symmetric(tau)
        n = (L - self.sr_of(tau)[1]) // self.nbb
        nodes = state.atoms
        on = set(path)
        rad = [a for a in atoms if nodes[a].element != "H" and nodes[a].radical]
        pis = [i for i in range(L - 1) if nodes[path[i]].edges.get(path[i + 1]) == 2.0]
        if len(rad) + len(pis) > 1 or any(a not in on for a in rad) or any(nodes[a].radical > 1 for a in rad):
            raise _OffDeclared("a product with more than one feature, or a feature off the backbone")
        x = side = None
        if rad:
            i = path.index(rad[0])
            if i in (0, L - 1):
                family, feature = "radical_end", (rad[0],)
                side = "start" if i == 0 and L > 1 and not sym else None
            else:
                family, feature, x = "radical_interior", (rad[0],), (min(i, L - 1 - i) if sym else i)
        elif pis:
            i = pis[0]
            if i == 0 or i + 1 == L - 1:
                family = "alkene_terminal"
                if self.legacy or i == 0 and sym:
                    feature = (path[0], path[1]) if i == 0 else (path[L - 1], path[L - 2])
                elif i + 1 == L - 1:
                    feature = (path[L - 1], path[L - 2])
                else:
                    feature, side = (path[0], path[1]), "start"
            else:
                family = "alkene_internal"
                if sym and not self.legacy and i > L - 2 - i:
                    feature, x = (path[i + 1], path[i]), L - 2 - i
                else:
                    feature, x = (path[i], path[i + 1]), (min(i, L - 2 - i) if sym else i)
        else:
            family, feature = "chain", ()
        label = self.make_label(family, state, feature, tau)
        at = x if family in INTERIOR_FAMILIES else side
        signature = (family, n, at, tau, me._comp_key(state, cid),
                     min(tuple(_label(state, [a]) for a in path), tuple(_label(state, [a]) for a in reversed(path))))
        if signature not in self._verified:
            key, forms = self.forms(family, n, at, tau)
            g = me._graph(me._first_rows(state, cid))
            if me._comp_key(state, cid) != key or not any(_isomorphic(g, h) for h in forms):
                raise _OffDeclared(f"a product is not the declared {family} of {n} periods")
            if family not in INTERIOR_FAMILIES and self.label(family, n, side, tau) != label:
                raise _OffDeclared(f"a product labelled {label}, the probe {self.label(family, n, side, tau)}")
            self._verified.add(signature)
        ends = tuple((seg[0] in atoms, seg[-1] in atoms, any(a in atoms for a in seg)) for seg in segs)
        place = None
        if family in INTERIOR_FAMILIES:
            place = _place(feature, segs, atoms, x)._replace(pos=_product_pos(x, L, WIDTH[family], sym),
                                                             flipped=flipped, oriented=not sym)
        return family, label, n, ends, place, L, (family if self.legacy else label)

    def identify(self, state, cid) -> Tuple:
        """(label, n, positions) of a component with no backbone origin (a
        depository molecule): the probe it is isomorphic to, its chain type
        and length from its molar mass and the family's mass offset per type."""
        mass = Fraction(state.components.molar_mass_g_mol(cid))
        g = me._graph(me._first_rows(state, cid))
        key = me._comp_key(state, cid)
        for family in FAMILIES:
            for tau in [None] + sorted(t for t in self.sr if t != self.declared):
                offs = self.offsets.get(family, ()) if tau is None else self.type_offsets.get((family, tau), ())
                r = self.sr_of(tau)[1]
                for off in sorted(offs):
                    n = (mass - off) / self.unit_mass
                    if n.denominator != 1 or n < (1 if r == 0 else 0):
                        continue
                    n = int(n)
                    hits = []
                    for x in self.positions(family, n, tau):
                        try:
                            k, forms = self.forms(family, n, x, tau)
                        except _OffDeclared:
                            continue
                        if k == key and any(_isomorphic(g, h) for h in forms):
                            hits.append(x)
                    if hits:
                        if tau is None and self.legacy:
                            label = self.label(family, n)
                        elif family in INTERIOR_FAMILIES:
                            label = self._label_at(family, n, hits[0], tau)
                        else:
                            label = self.label(family, n, hits[0], tau)
                        return label, n, (None if hits == [None] else hits)
        raise _OffDeclared("a depository molecule that is no declared state")


def _fired(engine, state, tr, view, segs, probes) -> List[Tuple]:
    """The classified product components of firing ``tr`` on ``state``
    (rolled back)."""
    saved, engine.state = engine.state, state
    try:
        tx, cids = engine._begin_apply(tr, view)
        try:
            return [probes.classify(state, c, segs) for c in cids]
        finally:
            tx.rollback()
    finally:
        engine.state = saved


def _relative(product, L, k, nbb) -> Tuple:
    """(token, alpha, beta, gamma) of a product of one reactant of L
    backbone atoms with anchor k: its length is alpha*L + beta*k + gamma
    backbone atoms (the piece holding the reactant's first atom, its last,
    both or neither); the token is its family (``kmc-kernel-rates/3`` state
    model) or its class label."""
    family, _label, n, ends = product[:4]
    first, last, _any = ends[0]
    lp, family = (product[5], product[6]) if len(product) > 5 else (n * nbb, family)
    if first and last:
        return (family, 1, 0, lp - L)
    if first:
        return (family, 0, 1, lp - k)
    if last:
        return (family, 1, -1, lp - L + k)
    return (family, 0, 0, lp)


def _absolute(product) -> Tuple:
    family, label, n, _ends = product[:4]
    return (label, n)


# --------------------------------------------------------------------------
# positions (I-272): where an interior product's feature sits
# --------------------------------------------------------------------------
WIDTH = {"radical_interior": 1, "alkene_internal": 2}     # backbone atoms a feature of the family spans


class _Place(NamedTuple):
    """Where an interior product's feature sits on the reactant backbone it
    came from: the index of its first atom on that segment (None when the
    feature is on no segment), the first and last index of the segment the
    product holds and whether those are the segment's ends, and the
    feature's end distance on the product; ``pos``, the feature's position
    code on the product (``_product_pos``), ``flipped``, whether the
    product's head-to-tail reading runs against the segment's, and
    ``oriented``, whether the product's chain type is not its own reverse."""
    index: Optional[int]
    lo: Optional[int]
    hi: Optional[int]
    first: bool
    last: bool
    width: int
    x: int
    pos: Optional[Tuple[str, int]] = None
    flipped: bool = False
    oriented: bool = False


def _cut_end(p: _Place, at_last: bool) -> str:
    """The position kind of a distance counted from a fragment's cut end, the segment's last held atom
    (``at_last``) or its first: ``end_distance`` on a chain type that is its own reverse, else the end of the
    product's head-to-tail reading that the cut end is."""
    if not p.oriented:
        return "end_distance"
    return "tail_distance" if at_last != p.flipped else "head_distance"


def _place(feature, segs, atoms, x) -> _Place:
    for seg in segs:
        at = [i for i, a in enumerate(seg) if a in feature]
        if len(at) == len(feature):
            held = [i for i, a in enumerate(seg) if a in atoms]
            return _Place(at[0], held[0], held[-1], held[0] == 0, held[-1] == len(seg) - 1, len(feature), x)
    return _Place(None, None, None, False, False, len(feature), x)


def _placed(product, L, k, nbb, fixed: bool) -> Tuple:
    """The positioned code of one product of a reactant of L backbone atoms
    at anchor k: its length code (``_relative``; a join as ``("join",
    family, gamma)``) and, for an interior feature, where it sits.
    ``fixed``: the anchor is at a fixed end distance (an end, an edge site,
    an end-distance law), so the position is the feature's own end distance
    on the product, ``("end_distance", d)``.  Otherwise (a bulk site; an
    interior reactant in the far field) the position moves with the anchor:
    ``("offset", c)``, the feature c backbone atoms past the anchor on a
    product holding the whole reactant, or ``("end_distance", d)`` from the
    cut end of a fragment holding one reactant end.  A join's feature is
    ``("junction", True)``: where it sits depends on both reactants.  On a
    product whose chain type is not its own reverse an end distance names
    its end (``head_distance``, ``tail_distance``: ``_product_pos``,
    ``_cut_end``)."""
    if product[0] == "join":
        return tuple(product[:3]) + ((("junction", True),) if product[3] is not None else ())
    base = _relative(product, L, k, nbb)
    p = product[4]
    if p is None:
        return base
    if fixed:
        return base + (p.pos,)
    if p.index is None:
        raise _OffDeclared("an interior product whose feature is on no reactant backbone")
    if p.first and p.last:
        if p.oriented and p.flipped:
            raise _OffDeclared("an interior product holding its whole reactant, read against its reactant")
        return base + (("offset", p.index - k),)
    if p.first:
        return base + ((_cut_end(p, True), p.hi - p.index - (p.width - 1)),)
    if p.last:
        return base + ((_cut_end(p, False), p.index - p.lo),)
    raise _OffDeclared("an interior product that holds neither end of its reactant")


def _place_of(product) -> Optional[_Place]:
    return product[3] if product[0] == "join" else product[4]


def _across_mirror(L, rows) -> bool:
    """Whether a measured interior state of L backbone atoms has a product
    read across the mirror of its piece: a product holding the whole
    reactant, on a chain type that is its own reverse, with its feature past
    the middle (nearer the last atom), so read from the other end than the
    reactant's feature (whose positions are those nearer the head,
    ``_Probes.positions``)."""
    for _k, prods in rows:
        for r in prods:
            p = None if r[0] == "join" else r[4]
            if p is not None and p.index is not None and p.first and p.last and not p.oriented \
                    and 2 * p.index > L - p.width:
                return True
    return False


def _absolute_placed(product) -> Tuple:
    """``_absolute`` with an interior feature's end distance on the product."""
    if product[0] == "join":
        return _placed(product, 0, 0, 0, True)
    base = _absolute(product)
    return base if product[4] is None else base + (product[4].pos or ("end_distance", product[4].x),)


def _placed_json(code) -> list:
    if code and isinstance(code[-1], tuple):
        kind, *v = code[-1]
        return list(code[:-1]) + [{kind: v[0] if len(v) == 1 else list(v)}]
    return list(code)


def _placed_list(counter: Counter) -> list:
    return [[m, [_placed_json(p) for p in outcome]] for outcome, m in sorted(counter.items(), key=repr)]


def _placed_counter(rows, place) -> Counter:
    """Counter of the sorted positioned products of measured rows [(k, products)], each product coded by
    ``place(product, k)``."""
    return Counter(tuple(sorted((place(r, k) for r in prods), key=repr)) for k, prods in rows)


def _bulk_place(code, lo, hi=None) -> Tuple:
    """A bulk-row code with its ``("offset", c)`` as the end-distance range
    ``("bulk", p, q)``: over the bulk anchors [lo, L - 1 - hi] (those of one
    phase: ``_phase_bounds``) the feature's first atom is uniform on the
    product's backbone atoms [p, Lp - w - q] that carry its class (Lp backbone
    atoms, the feature spanning w)."""
    hi = lo if hi is None else hi
    if not code or not isinstance(code[-1], tuple) or code[-1][0] != "offset":
        return code
    c = code[-1][1]
    return code[:-1] + (("bulk", lo + c, hi - c - (WIDTH[family_of(code[0])] - 1)),)


def _phase_bounds(edge, L, nbb, ph) -> Tuple[int, int]:
    """(the first bulk anchor, the atoms past the last) of phase ph (the anchors k with k % nbb == ph) over the
    bulk anchors [edge, L - 1 - edge] of a reactant of L backbone atoms."""
    first = edge + (ph - edge) % nbb
    last = (L - 1 - edge) - (L - 1 - edge - ph) % nbb
    return first, L - 1 - last


def _end_distance(d, lp, family) -> int:
    """The end distance on a product of lp backbone atoms of a feature d from one of its ends."""
    return min(d, lp - WIDTH[family] - d)


def _expand_placed(outcomes, L, k, nbb, types=None) -> Counter:
    """Absolute positioned outcomes of relative positioned ones (JSON) at
    reactant length L (backbone atoms) and anchor k: (family, n), (family,
    n, end distance) or (``join``, family, gamma).  With ``types`` (the
    table's ``chain_types``) a token is a class label and an interior
    product (label, n, position kind, d) (``_product_pos``)."""
    out = Counter()
    for m, prods in outcomes:
        coded = []
        for p in prods:
            p, where = (p[:-1], p[-1]) if p and isinstance(p[-1], dict) else (p, None)
            if p[0] == "join":
                coded.append(tuple(p))
                continue
            token, alpha, beta, gamma = p
            lp = alpha * L + beta * k + gamma
            r = _residue(token, types)
            if lp < 1 or lp < r or (lp - r) % nbb:
                raise _OffDeclared(f"a positioned law gives an impossible product at {L} backbone atoms")
            if where is None:
                coded.append((token, (lp - r) // nbb))
                continue
            (kind, v), = where.items()
            if kind == "bulk":
                raise KernelRatesError("a bulk position is expanded anchor by anchor (_predict_placed_far)")
            family = family_of(token)
            if types is None:
                d = v if kind == "end_distance" else k + v
                coded.append((family, lp // nbb, _end_distance(d, lp, family)))
            elif kind == "end_distance":
                coded.append((token, (lp - r) // nbb, kind, _end_distance(v, lp, family)))
            else:
                w = WIDTH[family]
                x = {"head_distance": v, "tail_distance": lp - w - v, "offset": k + v}[kind]
                coded.append((token, (lp - r) // nbb) + _product_pos(x, lp, w, _symmetric(token, types)))
        out[tuple(sorted(coded, key=repr))] += m
    return out


def _exact_placed(outcomes, types=None) -> Counter:
    """Positioned exact outcomes (JSON: labels, lengths, end distances) as
    (family, n[, end distance]) / (``join``, family, gamma) tuples; with
    ``types``, (label, n[, position kind, d])."""
    out = Counter()
    for m, prods in outcomes:
        coded = []
        for p in prods:
            p, where = (p[:-1], p[-1]) if p and isinstance(p[-1], dict) else (p, None)
            if p[0] == "join":
                coded.append(tuple(p))
            elif types is not None:
                coded.append((p[0], int(p[1])) + (() if where is None else next(iter(where.items()))))
            elif where is None or "end_distance" not in where:
                coded.append((family_of(p[0]), int(p[1])))
            else:
                coded.append((family_of(p[0]), int(p[1]), where["end_distance"]))
        out[tuple(sorted(coded, key=repr))] += m
    return out


def _predict_placed_far(far, anchor, n, nbb, x=None, types=None, r=0, side=None) -> Counter:
    L = n * nbb + r
    if anchor == "end":
        return _expand_placed(far["outcomes"], L, 0 if side == "start" else L - 1, nbb, types)
    if anchor in ("interior", "end_distance"):
        return _expand_placed(far["outcomes"], L, x, nbb, types)
    edge, total = far["edge"], Counter()
    for k in range(L):
        if k < edge:
            total.update(_expand_placed(far["first"][str(k)], L, k, nbb, types))
        elif L - 1 - k < edge:
            total.update(_expand_placed(far["last"][str(L - 1 - k)], L, k, nbb, types))
        else:
            lo = _phase_bounds(edge, L, nbb, k % nbb)[0]
            bulk = far["bulk_phases"][str(k % nbb)] if "bulk_phases" in far else far["bulk"]
            rows = [[m, [p[:-1] + [{"offset": p[-1]["bulk"][0] - lo}]
                         if p and isinstance(p[-1], dict) and "bulk" in p[-1] else p for p in prods]]
                    for m, prods in bulk]
            total.update(_expand_placed(rows, L, k, nbb, types))
    return total


def _at_key(x, L, family, types, token) -> str:
    """The ``at`` key of an end-distance law for its reactant's feature at position x on L backbone atoms: the
    end distance on a chain type that is its own reverse, else ``<kind>=<d>`` (``_oriented_pos``)."""
    if _symmetric(token, types):
        return str(x)
    kind, d = _oriented_pos(x, L, WIDTH[family])
    return f"{kind}={d}"


def predict_placed(law: dict, n: int, nbb: int, x: Optional[int] = None, types: Optional[dict] = None) -> Counter:
    """The absolute positioned outcomes (Counter of product tuples: (family,
    n), (family, n, end distance) for an interior feature, (``join``,
    family, gamma)) of a law on a reactant of n periods: its ``positioned``
    outcomes (exact up to ``exact_to_n``, its far rule beyond), or an
    end-distance law's outcomes with the reactant's feature at end distance
    ``x`` (``exact`` below ``from_n``, ``at`` from it).  ``types``: the
    table's ``chain_types`` (``_expand_placed``); x is then the feature's
    position from the reactant's head (its end distance on a type that is its
    own reverse)."""
    L = n * nbb + _residue(law["class"], types)
    if law["anchor"] == "end_distance":
        if n < int(law["from_n"]):
            entry = law["exact"].get(str(n), {}).get(str(x))
            if entry is not None and "refused" in entry:
                raise KernelRatesError(f"{law['class']} at end distance {x} on {n} periods: {entry['refused']}")
            return Counter() if entry is None else _exact_placed(entry["outcomes"], types)
        entry = law["at"].get(_at_key(x, L, law["family"], types, law["class"]))
        return Counter() if entry is None else _expand_placed(entry["outcomes"], L, x, nbb, types)
    placed = law["positioned"]
    if law["anchor"] != "interior" and n <= int(law.get("exact_to_n") or 0):
        return _exact_placed(placed["exact"].get(str(n), []), types)
    return _predict_placed_far(placed["far"], law["anchor"], n, nbb, x, types, L - n * nbb, law.get("side"))


def _far_at(law: dict, n: int, nbb: int, x: int, types: Optional[dict] = None) -> Counter:
    """The absolute positioned outcomes of an interior law's far rule with its feature at backbone atom x on n
    periods (its positioned far rule, or its plain one when no product has an interior feature)."""
    placed = law.get("positioned")
    return _expand_placed((placed or law)["far"]["outcomes"], n * nbb + _residue(law["class"], types), x, nbb,
                          types)


def predict_near(law: dict, near: Optional[dict], n: int, nbb: int, x: int, types: Optional[dict] = None) -> Counter:
    """The absolute positioned outcomes (``predict_placed``'s tuples) of an
    interior participant of a law term with its feature at canonical end
    distance x on n periods (``types``: ``predict_placed``): its near-end law
    (``law_near_end``) where it speaks (``exact`` below ``from_n``, every
    position; ``at`` from it), the far law elsewhere."""
    if near is not None:
        if n < int(near["from_n"]):
            entry = near["exact"].get(str(n), {}).get(str(x))
            if entry is not None and "refused" in entry:
                raise KernelRatesError(f"{near['class']} at end distance {x} on {n} periods: {entry['refused']}")
            return Counter() if entry is None else _exact_placed(entry["outcomes"], types)
        L = n * nbb + _residue(near["class"], types)
        entry = near["at"].get(_at_key(x, L, near["family"], types, near["class"]))
        if entry is not None:
            return _expand_placed(entry["outcomes"], L, x, nbb, types)
    return _far_at(law, n, nbb, x, types)


def _law_view(outcomes: Counter, by_class: bool = False) -> Counter:
    """Absolute outcomes at the resolution of a far rule's relative codes: an
    interior product by its family, since its class is read at its position
    (which the codes do not keep: one length of an interior family can carry
    several classes); ``by_class`` (codes by class label): unchanged."""
    if by_class:
        return outcomes
    out = Counter()
    for prods, m in outcomes.items():
        out[tuple(sorted(((family_of(p[0]), p[1]) if p[0] != "join" and family_of(p[0]) in INTERIOR_FAMILIES
                          else p) for p in prods))] += m
    return out


def _outcome_list(counter: Counter) -> list:
    return [[m, [list(p) for p in outcome]] for outcome, m in sorted(counter.items())]


class _Laws:
    """The product law of every kernel term (module doc, product laws)."""

    def __init__(self, engine, probes, n_hi):
        self.engine, self.probes, self.n_hi = engine, probes, n_hi
        self.nbb = probes.nbb
        self.types = probes.types_doc()
        self._partners = {}

    def _class(self, label):
        """(family, chain type or None, residue) of a class label."""
        tau = self.probes._tau(type_of(label))
        return family_of(label), tau, self.probes.sr_of(tau)[1]

    # -- the states of one participant class ---------------------------------
    def _states(self, label):
        """[(n, x)] the law reads for a participant class, every one a probe
        state carrying ``label``: every length of a chain or end family (an
        end family at the end carrying it), and the far-field interior
        positions (the two middlemost of the longest probe, the middlemost of
        the next)."""
        family, tau, r = self._class(label)
        probes = self.probes
        if family in INTERIOR_FAMILIES:     # the middle of the two longest probes
            longest, shorter = self._middle(label)
            return longest[:2] + shorter[:1]
        return [(n, x) for n in range(0 if r else 1, self.n_hi + 1) for x in probes.positions(family, n, tau)
                if probes.carries(label, n, x)]

    def _middle(self, label):
        """([(n, x)] on the longest probe, the same on the next): every
        interior state carrying ``label``, middlemost first."""
        family, tau, r = self._class(label)
        probes = self.probes
        if probes.symmetric(tau):
            order = lambda n: list(reversed(probes.positions(family, n, tau)))
        else:
            order = lambda n: sorted(probes.positions(family, n, tau),
                                     key=lambda x: (abs(2 * x + WIDTH[family] - n * self.nbb - r), x))
        mine = lambda n: [(n, x) for x in order(n) if probes.carries(label, n, x)]
        return mine(self.n_hi), mine(self.n_hi - 1)

    # -- measuring --------------------------------------------------------------
    def _uni(self, ci, label, n, x):
        """[(anchor, products)] of every allowed transition of class ``ci`` at
        the feature of probe (label, n, x)."""
        family, tau, _r = self._class(label)
        seg, feature, entry = self.probes.entry(family, n, x, tau=tau)
        allowed = _uni_all(entry).get(ci, ((), 0))[0]
        if not allowed:
            return seg, []
        state, paths = self.probes.melt()
        self.probes.put(state, paths[0], family, n, x, tau)
        cid = state.components.component_of(seg[0])
        view = me._LocalView(state, [cid])
        out = []
        for tr in allowed:
            touched = _heavy_touched(tr)
            if feature and not touched & set(feature):
                continue
            k = min(seg.index(a) for a in touched if a in seg)
            out.append((k, _fired(self.engine, state, tr, view, [seg], self.probes)))
        return seg, out

    def _partner(self, ci, q, label, markers=None):
        """(family, n, x, sites, chain type) of a probe with sites of
        participant ``q`` of class ``ci`` in class ``label`` (of a deferral
        marker in ``markers``, when given): the far field first."""
        key = (ci, q, label, markers)
        if key not in self._partners:
            family, tau, _r = self._class(label)
            found = None
            for n, x in sorted(self._states(label), key=lambda s: -s[0]) + \
                    [(self.n_hi, xx) for xx in self.probes.positions(family, self.n_hi, tau)
                     if self.probes.carries(label, self.n_hi, xx)]:
                seg, feature, entry = self.probes.entry(family, n, x, chain=1, tau=tau)
                sites = [(kind, s) for (cj, p, kind), ss in sorted(entry.sites.items(), key=repr)
                         if cj == ci and p == q and (markers is None or _marker(kind) in markers)
                         for s in ss if not feature or set(s) & set(feature)]
                if sites:
                    found = (family, n, x, sites, tau)
                    break
            if found is None:
                raise _OffDeclared(f"no partner site of class {label}")
            self._partners[key] = found
        return self._partners[key]

    def _bi(self, ci, p, label, n, x, partner_label, check=False, markers=None, partner_markers=None):
        """[(anchor, product)] of every site of participant ``p`` at the
        feature of probe (label, n, x), each fired with a partner of class
        ``partner_label``; the product is participant p's (``("join", token,
        gamma)`` when one component holds both)."""
        family, tau, _r = self._class(label)
        seg, feature, entry = self.probes.entry(family, n, x, tau=tau)
        cls = self.engine.classes[ci]
        mine = [(kind, s) for (cj, pp, kind), ss in sorted(entry.sites.items(), key=repr)
                if cj == ci and pp == p and (markers is None or _marker(kind) in markers)
                for s in ss if not feature or set(s) & set(feature)]
        if not mine:
            return seg, []
        q = 1 - p
        pf, pn, px, theirs, ptau = self._partner(ci, q, partner_label, partner_markers)
        state, paths = self.probes.melt()
        self.probes.put(state, paths[0], family, n, x, tau)
        pseg, _pfeat = self.probes.put(state, paths[1], pf, pn, px, ptau)
        ca, cb = state.components.component_of(seg[0]), state.components.component_of(pseg[0])
        out = []
        for kind, site in mine:
            usable = [(kb, sb) for kb, sb in theirs
                      if not me._pair_leak(cls, *((kind, kb) if p == 0 else (kb, kind)))]
            if not usable:
                raise _OffDeclared("a site whose every partner leaks")
            seen, kept = set(), None
            for kb, sb in usable[:2 if check else 1]:
                c0, s0, c1, s1 = (ca, site, cb, sb) if p == 0 else (cb, sb, ca, site)
                saved, self.engine.state = self.engine.state, state
                try:
                    tr, view = self.engine._pair_transition(ci, c0, s0, c1, s1)
                finally:
                    self.engine.state = saved
                if tr.leak:
                    raise _OffDeclared("a leaking pair")
                prods = _fired(self.engine, state, tr, view, [seg, pseg], self.probes)
                own = [r for r in prods if r[3][0][2]]
                if len(own) != 1:
                    raise _OffDeclared("a participant split over several products")
                (r,) = own
                if r[3][1][2]:
                    if not (r[3][0][0] and r[3][0][1] and r[3][1][0] and r[3][1][1]):
                        raise _OffDeclared("a joined product that is not both whole reactants")
                    got = ("join", r[6], r[5] - len(seg) - len(pseg))
                    full = got + (r[4],)
                else:
                    got, full = r[:4], r
                seen.add(got)
                kept = kept or full
            if len(seen) != 1:
                raise _OffDeclared("a participant's product depends on its partner's site")
            k = min(seg.index(a) for a in site if a in seg)
            out.append((k, kept))
        return seg, out

    # -- the law of one participant ------------------------------------------
    def participant(self, ci, p, label, partner_label=None, markers=None, partner_markers=None) -> dict:
        family, _tau, r_ = self._class(label)
        anchor = "chain" if family == "chain" else "end" if family in END_FAMILIES else "interior"
        nbb, by_class = self.nbb, self.types is not None
        raw = {}

        def measure(n, x, check=False):
            if partner_label is None:
                seg, rows = self._uni(ci, label, n, x)
                rows = [(k, tuple(sorted(prods, key=repr))) for k, prods in rows]
            else:
                seg, rows = self._bi(ci, p, label, n, x, partner_label, check, markers, partner_markers)
                rows = [(k, (prod,)) for k, prod in rows]
            L = len(seg)
            if anchor == "end":
                rows = [(0 if x == "start" else L - 1, r) for _k, r in rows]
            elif anchor == "interior":
                rows = [(x, r) for _k, r in rows]

            def code(prods, k, rel):
                out = []
                for r in prods:
                    if r[0] == "join":
                        out.append(r[:3])
                    else:
                        out.append(_relative(r, L, k, nbb) if rel else _absolute(r))
                return tuple(sorted(out))
            rel = [(k, code(prods, k, True)) for k, prods in rows]
            ab = Counter(code(prods, k, False) for k, prods in rows)
            raw[(n, x)] = (L, rows)
            return L, rel, ab

        # one length rule for every participant: the law reads the probe states
        # carrying its class, and its far rule needs the class at the law's reach
        # (and, along a chain or at an end, on EXACT_POINTS lengths to check it on)
        states = self._states(label)
        lengths = sorted({n for n, _x in states})
        # a class of an end family carried at short lengths only (a length label of the family, e.g. the
        # methyl radical) has no longer state: its law is exact on those lengths, with no far outcome
        fixed = anchor == "end" and bool(lengths) and lengths[-1] < self.n_hi
        if not lengths or (not fixed and (lengths[-1] != self.n_hi or
                                          (anchor != "interior" and len(lengths) < EXACT_POINTS))):
            raise _OffDeclared(f"no probe of class {label} at the law's reach ({self.n_hi} periods; "
                               f"probe lengths carrying it: {lengths or 'none'})")
        law = {"class": label, "family": family, "anchor": anchor}
        if anchor == "end" and any(x == "start" for _n, x in states):
            law["side"] = "start"
        if fixed:
            exact = {n: measure(n, x, check=True)[2] for n, x in states}
            law.update({"fixed_lengths": lengths, "far": {"outcomes": []},
                        "exact": {str(n): _outcome_list(c) for n, c in sorted(exact.items())},
                        "exact_to_n": lengths[-1]})
            return self._positioned(law, raw)
        if anchor == "interior":
            # the far field: the two middlemost states of the longest probe and the middlemost of the next,
            # each read off the mirror of a piece that is its own reverse (module doc, chain types)
            got = []
            for want, candidates in zip((2, 1), self._middle(label)):
                kept = 0
                for n, x in candidates:
                    if kept == want:
                        break
                    rel = measure(n, x, check=not got)[1]
                    if _across_mirror(*raw[(n, x)]):
                        del raw[(n, x)]
                        continue
                    got.append(Counter(r for _k, r in rel))
                    kept += 1
            if not got:
                raise _OffDeclared(f"{label}: no far-field state off the mirror of its piece")
            if any(g != got[0] for g in got[1:]):
                raise _OffDeclared(f"{label}: the far-field products of class {ci} differ between "
                                   "neighbouring positions")
            law["far"] = {"outcomes": _outcome_list(got[0])}
            return self._positioned(law, raw)
        exact, rel_hi = {}, {}
        for n, x in states:
            L, rel, ab = measure(n, x, check=(n == self.n_hi))
            exact[n] = ab
            rel_hi[n] = (L, rel)
        if anchor == "end":
            L, rel = rel_hi[self.n_hi]
            far = Counter(r for _k, r in rel)
            law["far"] = {"outcomes": _outcome_list(far)}
        else:
            L, rel = rel_hi[self.n_hi]
            by_k = defaultdict(Counter)
            for k, r in rel:
                by_k[k][r] += 1
            # the bulk is periodic in the anchor: one row per phase k % nbb, read at the middlemost anchor
            mid = L // 2
            bulk = {ph: by_k.get(mid + (ph - mid) % nbb, Counter()) for ph in range(nbb)}
            odd = [k for k in range(L) if by_k.get(k, Counter()) != bulk[k % nbb]]
            edge = max([min(k, L - 1 - k) + 1 for k in odd], default=0)
            if 2 * edge >= L:
                raise _OffDeclared(f"class {ci}: end effects reach the middle of the longest probe")
            rows = {ph: _outcome_list(c) for ph, c in bulk.items()}
            law["far"] = {"edge": edge,
                          "first": {str(k): _outcome_list(by_k.get(k, Counter())) for k in range(edge)},
                          "last": {str(k): _outcome_list(by_k.get(L - 1 - k, Counter())) for k in range(edge)}}
            if all(v == rows[0] for v in rows.values()):
                law["far"]["bulk"] = rows[0]
            else:
                law["far"]["bulk_phases"] = {str(ph): v for ph, v in rows.items()}
        for n in sorted(exact)[-EXACT_POINTS:]:
            if len(exact) < EXACT_POINTS or \
                    _law_view(predict(law, n, nbb, self.probes.label, self.types), by_class) != \
                    _law_view(exact[n], by_class):
                raise _OffDeclared(f"class {ci} {label}: the far-field law does not hold at {n} periods")
        law["exact"] = {str(n): _outcome_list(c) for n, c in sorted(exact.items())}
        law["exact_to_n"] = max(exact) if exact else 0
        return self._positioned(law, raw)

    # -- positions (I-272) ------------------------------------------------------
    def _positioned(self, law, raw) -> dict:
        """``law`` with ``positioned`` (module doc, positions) when any product
        has an interior feature: the law's outcomes resolved by where each
        feature sits, read from the same fired transitions; its far rule must
        reproduce the positioned exact outcomes on the trailing lengths, or
        the law says ``positions_refused`` and why."""
        if not any(_place_of(r) is not None for _L, rows in raw.values() for _k, prods in rows for r in prods):
            return law
        try:
            law["positioned"] = self._placed_law(law, raw)
        except _OffDeclared as error:
            law["positions_refused"] = str(error)
        return law

    def _placed_law(self, law, raw) -> dict:
        nbb, anchor = self.nbb, law["anchor"]
        out = {}
        if anchor == "interior":
            got = [_placed_counter(rows, lambda r, k, L=L: _placed(r, L, k, nbb, False))
                   for L, rows in raw.values()]
            if any(g != got[0] for g in got[1:]):
                raise _OffDeclared("the far-field positions differ between neighbouring positions")
            return {"far": {"outcomes": _placed_list(got[0])}}
        out["exact"] = {str(n): _placed_list(_placed_counter(rows, lambda r, k: _absolute_placed(r)))
                        for (n, _x), (_L, rows) in sorted(raw.items())}
        if "fixed_lengths" in law:
            out["far"] = {"outcomes": []}
            return out
        L, rows = next(v for (n, _x), v in raw.items() if n == self.n_hi)
        if anchor == "end":
            out["far"] = {"outcomes": _placed_list(_placed_counter(rows, lambda r, k: _placed(r, L, k, nbb, True)))}
        else:
            by_k = defaultdict(list)
            for k, prods in rows:
                by_k[k].append((k, prods))
            moving = {k: _placed_counter(rs, lambda r, kk: _placed(r, L, kk, nbb, False)) for k, rs in by_k.items()}
            mid = L // 2
            bulk = {ph: moving.get(mid + (ph - mid) % nbb, Counter()) for ph in range(nbb)}
            odd = [k for k in range(L) if moving.get(k, Counter()) != bulk[k % nbb]]
            edge = max([min(k, L - 1 - k) + 1 for k in odd], default=0)
            if 2 * edge >= L:
                raise _OffDeclared("end effects on the positions reach the middle of the longest probe")
            fixed = lambda k: _placed_list(_placed_counter(by_k.get(k, []), lambda r, kk: _placed(r, L, kk, nbb, True)))
            placed = {}
            for ph, moved in bulk.items():
                lo, hi = _phase_bounds(edge, L, nbb, ph)
                placed[ph] = _placed_list(Counter({tuple(_bulk_place(p, lo, hi) for p in outcome): m
                                                   for outcome, m in moved.items()}))
            out["far"] = {"edge": edge,
                          "first": {str(k): fixed(k) for k in range(edge)},
                          "last": {str(k): fixed(L - 1 - k) for k in range(edge)}}
            if all(v == placed[0] for v in placed.values()):
                out["far"]["bulk"] = placed[0]
            else:
                out["far"]["bulk_phases"] = {str(ph): v for ph, v in placed.items()}
        r_ = self._class(law["class"])[2]
        for n in sorted(int(n) for n in law["exact"])[-EXACT_POINTS:]:
            if _predict_placed_far(out["far"], anchor, n, nbb, None, self.types, r_, law.get("side")) != \
                    _exact_placed(out["exact"][str(n)], self.types):
                raise _OffDeclared(f"the positioned far rule does not hold at {n} periods")
        return out

    def _position(self, pos, family, L) -> Tuple[int, str]:
        """(backbone atom from the head, ``at`` key) of a profile position: an end distance on a chain type
        that is its own reverse, else (kind, d) from the head or the tail."""
        if isinstance(pos, int):
            return pos, str(pos)
        kind, d = pos
        return (d if kind == "head_distance" else L - WIDTH[family] - d), f"{kind}={d}"

    # -- the law of an interior participant with sites near an end only -----------
    def participant_at(self, ci, p, label, profile, partner_label=None, weight=None, markers=None,
                       partner_markers=None) -> dict:
        """The law of an interior participant whose sites lie only near a
        chain end (``per_feature`` 0; module doc, end-distance laws): per end
        distance d with sites on chains of ``from_n`` periods and more, its
        sites and positioned outcomes on the longest probe with the feature at
        d (the count checked against the participant's ``profile``); on every
        shorter length, the same at every position (``exact``).  ``markers``,
        ``partner_markers``: those of the sites a deferring term counts (``_bi``)."""
        family, tau, r_ = self._class(label)
        nbb = self.nbb
        far, from_n = self._profile_reach(profile, nonzero=True, tau=tau)

        def rows_at(n, x, check=False):
            return self._rows_at(ci, p, label, n, x, partner_label, check, markers, partner_markers)

        at = {}
        for pos in sorted(far):
            x, key = self._position(pos, family, self.n_hi * nbb + r_)
            if not self.probes.carries(label, self.n_hi, x):
                continue
            L, rows = rows_at(self.n_hi, x, check=True)
            if not rows:
                continue
            if len(rows) not in far[pos]:
                raise _OffDeclared(f"{label} at end distance {key}: {len(rows)} sites on the probe, its profile "
                                   f"counts {far[pos]}")
            at[key] = {"sites": len(rows),
                       "outcomes": _placed_list(_placed_counter(rows, lambda r, k: _placed(r, L, x, nbb, True)))}
            if weight is not None:
                at[key]["k_per_feature"] = weight * len(rows)
        return {"class": label, "family": family, "anchor": "end_distance", "from_n": from_n, "at": at,
                "exact": self._exact_positions(family, label, from_n, rows_at, weight)}

    def _rows_at(self, ci, p, label, n, x, partner_label=None, check=False, markers=None, partner_markers=None):
        """(backbone atoms, [(anchor, products)]) of participant p's sites with its feature at x on n periods."""
        if partner_label is None:
            seg, rows = self._uni(ci, label, n, x)
            return len(seg), [(k, tuple(sorted(prods, key=repr))) for k, prods in rows]
        seg, rows = self._bi(ci, p, label, n, x, partner_label, check, markers, partner_markers)
        return len(seg), [(k, (prod,)) for k, prod in rows]

    def _profile_reach(self, profile, nonzero=False, tau=None):
        """({position: counts} of an interior participant's profile on the longest probe, from_n): from_n the
        shortest length from which every longer probe's profile is the same; a position is an end distance on
        a chain type that is its own reverse, else (``head_distance`` | ``tail_distance``, d).  ``nonzero``:
        positions with some site only."""
        nbb, ppu = self.nbb, self.probes.ppu
        r = self.probes.sr_of(tau)[1]
        counted = defaultdict(dict)                  # backbone atoms -> {position: counts}
        for key, values in profile.items():
            m = re.fullmatch(r"chain_backbone_atoms=(\d+),(end|head|tail)_distance=(\d+)", key)
            if m and (any(values) or not nonzero):
                pos = int(m[3]) if m[2] == "end" else (f"{m[2]}_distance", int(m[3]))
                counted[int(m[1])][pos] = sorted(values)
        lengths = [u * ppu * nbb + r for u in range(1, self.n_hi // ppu + 1)]
        far = counted.get(lengths[-1], {})
        from_nb = lengths[-1]
        for nb in reversed(lengths):
            if counted.get(nb, {}) != far:
                break
            from_nb = nb
        return far, -(-(from_nb - r) // nbb)

    def _exact_positions(self, family, label, from_n, rows_at, weight) -> dict:
        """Per length shorter than from_n, the sites and positioned outcomes at every canonical position with
        sites (an entry the probes cannot read says ``refused``)."""
        _family, tau, r = self._class(label)
        exact = {}
        for n in range(0 if r else 1, from_n):
            for x in self.probes.positions(family, n, tau):
                if not self.probes.carries(label, n, x):
                    continue
                try:
                    _L, rows = rows_at(n, x)
                except _OffDeclared as error:
                    exact.setdefault(str(n), {})[str(x)] = {"refused": str(error)}
                    continue
                if rows:
                    entry = {"sites": len(rows),
                             "outcomes": _placed_list(_placed_counter(rows, lambda r, k: _absolute_placed(r)))}
                    if weight is not None:
                        entry["k_per_feature"] = weight * len(rows)
                    exact.setdefault(str(n), {})[str(x)] = entry
        return exact

    # -- the near-end law of an interior participant with far-field sites (I-273) -------------------------
    def participant_near(self, ci, p, label, count, law, partner_label=None, weight=None, markers=None,
                         partner_markers=None) -> Optional[dict]:
        """The near-end law of an interior participant of a law term (module
        doc, near-end laws), or None where none is needed: on the longest
        probe, every end distance d where the participant's sites or their
        positioned outcomes differ from its far law at d (``at``, zero sites
        included), scanned one past the reach of its profile there (an
        internal alkene's profile key is one past its mirror bond's end
        distance; on a chain type that is not its own reverse, from each end);
        on every length shorter than where the profile settles, the sites and
        outcomes at every canonical position (``exact``).
        ``markers``, ``partner_markers``: ``participant_at``."""
        family, tau, r_ = self._class(label)
        nbb, L_hi = self.nbb, self.n_hi * self.nbb + r_
        far_profile, from_n = self._profile_reach(count.get("profile", {}), tau=tau)

        def rows_at(n, x, check=False):
            return self._rows_at(ci, p, label, n, x, partner_label, check, markers, partner_markers)

        if self.probes.symmetric(tau):
            reach = max(far_profile, default=0) + 1
            scan = [(d, d, str(d), reach) for d in self.probes.positions(family, self.n_hi, tau) if d <= reach]
        else:
            reach = {kind: max([d for kk, d in far_profile if kk == kind], default=0) + 1
                     for kind in ("head_distance", "tail_distance")}
            scan = []
            for x in self.probes.positions(family, self.n_hi, tau):
                kind, d = _oriented_pos(x, L_hi, WIDTH[family])
                if d <= reach[kind]:
                    scan.append((x, d, f"{kind}={d}", reach[kind]))
        at = {}
        for x, d, key, reach_d in scan:
            if not self.probes.carries(label, self.n_hi, x):
                continue
            L, rows = rows_at(self.n_hi, x, check=True)
            got = _exact_placed(_placed_list(_placed_counter(rows, lambda r, k: _absolute_placed(r))), self.types)
            if len(rows) == count["per_feature"] and got == _far_at(law, self.n_hi, nbb, x, self.types):
                continue
            if d == reach_d:
                raise _OffDeclared(f"{label}: its near-end law differs from the far law at end distance {key}, past "
                                   f"its profile's reach on {L_hi} backbone atoms")
            at[key] = {"sites": len(rows),
                       "outcomes": _placed_list(_placed_counter(rows, lambda r, k: _placed(r, L, x, nbb, True)))}
            if weight is not None:
                at[key]["k_per_feature"] = weight * len(rows)
        if not at and from_n <= 1:
            return None
        return {"class": label, "family": family, "anchor": "end_distance", "from_n": from_n, "at": at,
                "exact": self._exact_positions(family, label, from_n, rows_at, weight)}

    def near_end_law(self, ci, combo, order, participants, law, weight, markers=None) -> Optional[dict]:
        """``participant_near`` of every interior participant of a law term:
        the law itself for one participant, else per slot (``S0``/``S1``)
        with a near-end law; None when no participant needs one.
        ``markers``: per participant those of the sites a deferring term counts."""
        def one(p, slot_law, partner):
            c = participants[p].get(combo[p])
            if family_of(combo[p]) not in INTERIOR_FAMILIES or not c or not c.get("per_feature"):
                return None
            return self.participant_near(ci, p, combo[p], c, slot_law, partner, weight if partner is None else None,
                                         *((markers[p], markers[1 - p]) if markers else ()))
        if len(combo) == 1:
            return one(0, law, None)
        out = {f"S{slot}": one(p, law[f"S{slot}"], combo[1 - p]) for slot, p in enumerate(order)}
        out = {slot: near for slot, near in out.items() if near is not None}
        return out or None

    def end_distance_law(self, ci, combo, order, participants, weight, markers=None) -> dict:
        """The law of a term whose interior participant has sites only near
        a chain end: ``participant_at`` for every such participant, the far-
        field law (``participant``, its partner read at the end distance it
        has sites at) for the other.  ``markers``: ``near_end_law``."""
        def one(p, partner):
            c = participants[p][combo[p]]
            mark = (markers[p], markers[1 - p]) if markers else ()
            if family_of(combo[p]) in INTERIOR_FAMILIES and not c.get("per_feature"):
                return self.participant_at(ci, p, combo[p], c["profile"], partner, weight if partner is None else None,
                                           *mark)
            return self.participant(ci, p, combo[p], partner, *mark)
        if len(combo) == 1:
            return one(0, None)
        return {f"S{slot}": one(p, combo[1 - p]) for slot, p in enumerate(order)}


def _expand(outcomes, L, k, nbb, label_of, types=None) -> Counter:
    """Absolute outcomes of relative ones at reactant length L (backbone
    atoms) and anchor k; joins pass through; ``types``: codes by class label
    (``predict``)."""
    out = Counter()
    for m, prods in outcomes:
        coded = []
        for p in prods:
            if p[0] == "join":
                coded.append(tuple(p))
                continue
            token, alpha, beta, gamma = p
            lp = alpha * L + beta * k + gamma
            if types is None:
                if lp < 1 or lp % nbb:
                    return None
                coded.append((label_of(token, lp // nbb), lp // nbb))
                continue
            r = _residue(token, types)
            if lp < 1 or lp < r or (lp - r) % nbb:
                return None
            coded.append((token, (lp - r) // nbb))
        out[tuple(sorted(coded))] += m
    return out


def predict(law: dict, n: int, nbb: int, label_of=None, types: Optional[dict] = None) -> Counter:
    """The absolute outcomes (Counter of product tuples) of a chain or end
    law's far rule on n periods.  ``types``: the table's ``chain_types``
    (codes by class label, a reactant of its chain type's residue, a bulk
    row per phase)."""
    label_of = label_of or (lambda family, m: family)
    L = n * nbb + _residue(law["class"], types)
    far = law["far"]
    total = Counter()
    if law["anchor"] == "end":
        got = _expand(far["outcomes"], L, 0 if law.get("side") == "start" else L - 1, nbb, label_of, types)
        return Counter() if got is None else got
    edge = far["edge"]
    for k in range(L):
        if k < edge:
            rows = far["first"][str(k)]
        elif L - 1 - k < edge:
            rows = far["last"][str(L - 1 - k)]
        else:
            rows = far["bulk_phases"][str(k % nbb)] if "bulk_phases" in far else far["bulk"]
        got = _expand(rows, L, k, nbb, label_of, types)
        if got is None:
            raise _OffDeclared(f"a far-field law gives an impossible product at {n} periods")
        total.update(got)
    return total


def _template_overrides(engine, probes, t, ci, entries) -> list:
    """The depository overrides of one template at their exact sizes:
    per override its weight, its reactant states (label, n, interior
    positions or None; unordered for a bimolecular term), the products of a
    matched transition and the matched multiplicity (transitions on the
    reactant; for two reactants, ordered cross-component site pairs, the
    engine's count before its divisor), and the catalog's depository entry
    (``entries``, in override order; None where the catalog names none).
    Read on the override's own molecules
    with the melt engine's matching (``_uni_overrides``, ``_pair_override``).
    An override whose molecules are no declared state is kept with its
    reason (``refused``)."""
    if len(entries) != len(t.overrides):
        raise KernelRatesError(f"template {t.name}: {len(t.overrides)} overrides loaded, "
                               f"{len(entries)} in its payload")
    out = []
    for o, ov in enumerate(t.overrides):
        try:
            out.append(dict(_one_override(engine, probes, t, ci, o, ov), entry=entries[o]))
        except _OffDeclared as error:
            out.append({"name": ov.name, "refused": str(error), "entry": entries[o]})
    return out


def _one_override(engine, probes, t, ci, o, ov) -> dict:
    """One override of class ``ci`` (``_template_overrides``)."""
    state = me._rows_state(ov.reactants)
    comps = list(state.components.ids())
    if len(comps) != len(t.participants):
        raise _OffDeclared(f"override {ov.name}: {len(comps)} reactants for a "
                           f"{len(t.participants)}-participant template")
    reactants = [probes.identify(state, c) for c in comps]
    entries = engine._scan(comps, state=state)
    saved, engine.state = engine.state, state
    try:
        if len(comps) == 1:
            allowed = entries[comps[0]].uni_o.get((ci, o), ((), 0))[0]
            m = len(allowed)
            fire = (allowed[0], me._LocalView(state, comps)) if allowed else None
        else:
            m, fire = 0, None
            for ca, cb in ((comps[0], comps[1]), (comps[1], comps[0])):
                for (cj, p0, k0), s0 in sorted(entries[ca].sites.items(), key=repr):
                    if cj != ci or p0 != 0:
                        continue
                    for (ck, p1, k1), s1 in sorted(entries[cb].sites.items(), key=repr):
                        if ck == ci and p1 == 1 and engine._pair_override(ci, k0, k1) == o \
                                and not me._suppressed(k0, k1):
                            m += len(s0) * len(s1)
                            if fire is None:
                                fire = engine._pair_transition(ci, ca, s0[0], cb, s1[0])
    finally:
        engine.state = saved
    if not m:
        raise _OffDeclared(f"override {ov.name}: no transition of the template is its reaction")
    tr, view = fire
    saved, engine.state = engine.state, state
    try:
        tx, cids = engine._begin_apply(tr, view)
        try:
            products = sorted([list(probes.identify(state, c)[:2]) for c in cids])
        finally:
            tx.rollback()
    finally:
        engine.state = saved
    return {"name": ov.name, "weight": float(ov.weight), "multiplicity": m,
            "reactants": sorted([lab, n, pos] for lab, n, pos in reactants),
            "products": products}


# --------------------------------------------------------------------------
# the table
# --------------------------------------------------------------------------
def _load(doc, spec, temperature, allow_unbound):
    from rmgpy.kmc import run as kr
    from rmgpy.polymer_spec import LinearHomopolymerSpec
    if not isinstance(spec, LinearHomopolymerSpec):
        raise KernelRatesError("kernel rates need a linear homopolymer declaration (copolymer moments are "
                               "not defined on one primitive period)")
    if doc.get("schema") != "kmc-template-catalog/1" or not isinstance(doc.get("payloads"), list):
        raise KernelRatesError("not a kmc-template-catalog/1 document")
    bound = (doc.get("provenance") or {}).get("chemistry_sha256")
    if bound != kr.chemistry_sha256(spec) and not allow_unbound:
        raise KernelRatesError(f"the catalog's chemistry digest {bound!r} is not the declaration's "
                               f"({kr.chemistry_sha256(spec)})")


def _progress(message: str) -> None:
    print(f"kernel_rates: {message}", file=sys.stderr, flush=True)


def _law_term(laws, engine, probes, t, term, ci, combo, order, measured, placed_overrides, entries,
              participants, markers=None) -> None:
    """The product law (or why none) and the placed overrides of one kernel term;
    ``markers``: per participant the deferral markers of the sites it counts
    (None: every site)."""
    counts = [term] if len(combo) == 1 else [term["S0"], term["S1"]]
    if markers:
        # a deferring term's sites are not its row's: its laws count the term's own (I-275)
        participants = [{combo[p]: counts[order.index(p)]} for p in range(len(combo))]
    if any(family_of(c) in INTERIOR_FAMILIES and (x.get("per_feature") == 0 or x.get("k_per_feature") == 0)
           for c, x in zip([combo[p] for p in order], counts)):
        # an interior participant with no far-field site: the kernel rate is 0; its sites near an end have a
        # law per end distance (I-272)
        term["law"], term["law_inert"] = None, "an interior participant with no far-field site"
        term["overrides"] = []
        try:
            term["law_end_distance"] = laws.end_distance_law(ci, combo, order, participants, t.weight, markers)
        except _OffDeclared as error:
            term["law_end_distance_refused"] = str(error)
        return
    try:
        if len(combo) == 1:
            term["law"] = laws.participant(ci, 0, combo[0])
        else:
            term["law"] = {f"S{slot}": laws.participant(ci, p, combo[p], combo[1 - p],
                                                        *((markers[p], markers[1 - p]) if markers else ()))
                           for slot, p in enumerate(order)}
    except _OffDeclared as error:
        term["law"], term["law_refused"] = None, str(error)
    if term["law"] is not None:
        # an interior participant's sites near an end: which of them are left, and their products (I-273)
        try:
            near = laws.near_end_law(ci, combo, order, participants, term["law"], t.weight, markers)
        except _OffDeclared as error:
            term["law_near_end_refused"] = str(error)
        else:
            if near is not None:
                term["law_near_end"] = near
    if ci not in measured:
        measured[ci] = _template_overrides(engine, probes, t, ci, entries[ci])
    term["overrides"] = [o for o in measured[ci] if "refused" not in o
                         and sorted(r[0] for r in o["reactants"]) == sorted(combo)]
    placed_overrides.update((t.name, o["name"]) for o in term["overrides"])


def _chain_types_doc(probes, n_hi: int) -> dict:
    """``chain_types`` (module doc, chain types): per type its start offset
    and residue, whether it is its own reverse and whether it is the declared
    type; per family on it the classes it carries (an end family per end:
    ``last``, and ``start`` on a type that is not its own reverse, its class
    far from the ends and the length labels that differ; an interior family:
    the classes of its positions on the longest probe), its shortest length,
    its molar-mass and composition offsets (M = n * unit mass + offset) and
    its declared species."""
    out = {}
    for name, (s, r) in sorted(probes.sr.items()):
        tau = probes._tau(name)
        families = {}
        for family in FAMILIES:
            entry = {}
            if family in INTERIOR_FAMILIES:
                lengths = []
                for n in range(0 if r else 1, n_hi + 1):
                    xs = probes.positions(family, n, tau)
                    try:
                        if xs:
                            probes.forms(family, n, xs[0], tau)
                            lengths.append(n)
                    except _OffDeclared:
                        continue
                classes = set()
                for x in probes.positions(family, n_hi, tau):
                    try:
                        classes.add(probes._label_at(family, n_hi, x, tau))
                    except _OffDeclared:
                        continue
                if not lengths or not classes:
                    continue
                entry = {"classes": sorted(classes), "min_n": min(lengths)}
            else:
                ends = {}
                for x in probes.positions(family, n_hi, tau):
                    labels = {}
                    for n in range(0 if r else 1, n_hi + 1):
                        try:
                            labels[n] = probes.label(family, n, x, tau)
                        except _OffDeclared:
                            continue
                    if labels:
                        far = labels[max(labels)]
                        ends["start" if x == "start" else "last"] = {
                            "label_far": far, "min_n": min(labels),
                            "labels": {str(n): v for n, v in sorted(labels.items()) if v != far}}
                if not ends:
                    continue
                entry = {"ends": ends, "min_n": min(e["min_n"] for e in ends.values())}
            offsets = probes.offsets.get(family, set()) if tau is None else probes.type_offsets.get((family, tau), set())
            compositions = probes.composition_offsets.get(family, set()) if tau is None else \
                probes.type_composition_offsets.get((family, tau), set())
            if len(offsets) != 1 or len(compositions) != 1:
                raise KernelRatesError(f"{family} of chain type {name}: the molar mass or composition is not affine "
                                       f"in n ({sorted(offsets)}, {sorted(compositions)})")
            entry.update({"molar_mass_offset_g_mol": float(next(iter(offsets))),
                          "composition_offset": {el: str(v) for el, v in next(iter(compositions))},
                          "declared_species": {f"{lab}@{n}": sp for (lab, n), sp in sorted(probes.species.items())
                                               if family_of(lab) == family and type_of(lab) == tau}})
            families[family] = entry
        out[name] = {"start_offset": s, "residue": r, "symmetric": probes.sym[name],
                     "declared": name == probes.declared, "families": families}
    return out


def _unreached(kernels: dict, probes) -> dict:
    """The participant classes no term's law (``law``) reads, and the chain
    types none of whose classes a law reads."""
    participants = {c for k in kernels.values() for c in k["classes"]}
    lawful = {c for k in kernels.values() for t in k["templates"] if t["law"] is not None for c in k["classes"]}
    reached = {type_of(c) or probes.declared for c in lawful}
    return {"classes": sorted(participants - lawful), "chain_types": sorted(set(probes.sr) - reached)}


def kernel_rates(doc: dict, spec, temperature: float, *, run: Optional[dict] = None,
                 allow_unbound: bool = False, max_units: Optional[int] = None) -> dict:
    """The kernel-rate table (module doc) of a catalog document for one
    declaration at one temperature."""
    _load(doc, spec, temperature, allow_unbound)
    catalog = me.PayloadCatalog(doc["payloads"])
    geometry = ChainGeometry.of(spec)
    ppu = geometry.d_rep // geometry.period
    if ppu < 1 or geometry.d_rep % geometry.period:
        raise KernelRatesError(f"declared unit is not a whole number of primitive periods "
                               f"(d_rep {geometry.d_rep}, period {geometry.period})")
    templates = list(catalog.templates(float(temperature)))
    # per loaded template its overrides' depository entries (``loaded_templates`` keeps payload order)
    entries = [[o.get("entry") for o in t.get("overrides") or ()] for p in doc["payloads"]
               for t in p.get("templates", []) if t.get("status") == "emitted"]
    if len(entries) != len(templates):
        raise KernelRatesError(f"{len(templates)} templates loaded from {len(entries)} emitted ones")
    window = max(sum(1 for a in t.atoms if a.element != "H") for t in templates) if templates else 1
    bb_per_unit = len(geometry.backbone_keys)
    units_hi = max_units or max(2 * EXACT_POINTS + 2, math.ceil((2 * window + EXTRA_UNITS) / bb_per_unit) + 1)
    engine = me.MeltEngine(_single_chain_spec(spec, units_hi), catalog, float(temperature),
                           numpy.random.default_rng(0))
    if [t.name for t in engine.templates] != [t.name for t in templates]:
        raise KernelRatesError("the engine loaded another template set")
    counter = _Counter(engine)
    n_values = list(range(1, units_hi + 1))
    if bb_per_unit % ppu:
        raise KernelRatesError("the declared unit's backbone is not a whole number of primitive periods")
    nbb = bb_per_unit // ppu
    n_hi = units_hi * ppu
    probes = _Probes(spec, engine, ppu, nbb, 2 * units_hi + 3)
    chains, features = _states(spec, n_values, counter, ppu, nbb, probes)
    classes = engine.classes

    # per template, participant -> {class: count description}
    by_template = {ci: defaultdict(dict) for ci in range(len(templates))}
    kinds = {ci: defaultdict(lambda: defaultdict(set)) for ci in range(len(templates))}
    for name, chain in chains.items():
        for (ci, p, kind), series in chain.items():
            if any(series.values()):
                by_template[ci][p].setdefault(name, Counter())
                kinds[ci][p][name].add(kind)

    def chain_count(ci, p, markers=None, name="chain"):
        """The chain polynomial of participant ``p`` of class ``ci`` on the
        chain class ``name`` (its sites of a deferral marker in ``markers``,
        when given)."""
        total = Counter()
        for (cj, q, kind), series in chains[name].items():
            if cj == ci and q == p and (markers is None or _marker(kind) in markers):
                for n, v in series.items():
                    total[n] += v
        return _affine(dict(total))

    def feature_count(name, ci, p, markers=None):
        """The count of feature ``name`` (``chain_count``)."""
        mine = lambda k: k[0] == ci and k[1] == p and (markers is None or _marker(k[2]) in markers)
        c_far = sum(v for k, v in fars[name].items() if mine(k))
        prof = {}
        for position, counts in sorted(features[name]["profile"].items()):
            values = sorted({sum(v for k, v in c.items() if mine(k)) for c in counts})
            if values != [c_far]:
                prof[",".join(f"{position[i]}={position[i + 1]}" for i in range(0, len(position), 2))] = values
        return {"per_feature": c_far, "profile": prof}

    chain_poly = {}
    for ci in range(len(templates)):
        for p in list(by_template[ci]):
            for name in chains:
                if name in by_template[ci][p]:
                    chain_poly[(ci, p, name)] = chain_count(ci, p, name=name)
    feature_counts, fars = {}, {}
    for name, slot in features.items():
        far_by = slot["far"]
        if far_by:
            last = sorted(far_by)[-EXACT_POINTS:]
            if len(last) < EXACT_POINTS or any(far_by[n] != far_by[last[-1]] for n in last):
                raise KernelRatesError(f"feature {name}: the counts of a feature in the middle of the chain "
                                       "differ on the longest chains")
            far = far_by[last[-1]]
        else:
            far = Counter()       # a feature only near an end or on a short fragment: no far-field count
        fars[name] = far
        keys = set(far)
        for counts in slot["profile"].values():
            for c in counts:
                keys.update(c)
        for ci, p in sorted({(k[0], k[1]) for k in keys}):
            for k in keys:
                if k[0] == ci and k[1] == p:
                    kinds[ci][p][name].add(k[2])
            count = feature_count(name, ci, p)
            if count["per_feature"] or any(any(v) for v in count["profile"].values()):
                feature_counts[(ci, p, name)] = count
                by_template[ci][p][name] = Counter()

    kernels: Dict[str, dict] = {}
    kmc_only: List[dict] = []
    jobs = []                 # (term, class index, classes, role order, site markers): the product laws
    template_rows = []
    for ci, (t, cls) in enumerate(zip(templates, classes)):
        ktype = kernel_type(t, cls)
        row = {"template": t.name, "family": t.family, "type": ktype, "weight": float(t.weight),
               "weight_units": "m3/(mol*s) per matched transition" if cls.bimolecular
               else "1/s per matched transition",
               "divisor": cls.divisor, "phi_thinned": bool(cls.met_reverse), "participants": [],
               "kernels": []}
        parts = len(t.participants)
        for p in range(parts):
            entry = {}
            for name in sorted(by_template[ci].get(p, {})):
                if name in chains:
                    entry[name] = chain_poly[(ci, p, name)]
                else:
                    entry[name] = feature_counts[(ci, p, name)]
            row["participants"].append(entry)
        template_rows.append(row)
        if ktype is None:
            kmc_only.append({"template": t.name, "family": t.family, "reason": "no_kernel"})
            continue
        if any(not row["participants"][p] for p in range(parts)):
            empty = [p for p in range(parts) if not row["participants"][p]]
            if not cls.bimolecular and counter.leak[ci]:
                reason = "branching"
            elif any(not any(a.kind == "port" for a in t.atoms if a.participant == p) for p in empty):
                reason = "small_species_participant"     # a whole molecule that no declared chain piece is
            else:
                reason = "no_site_on_declared_states"
            kmc_only.append({"template": t.name, "family": t.family, "reason": reason})
            continue
        combos = [(c,) for c in row["participants"][0]] if parts == 1 else \
            [(c0, c1) for c0 in row["participants"][0] for c1 in row["participants"][1]]
        # per class pairing the deferral markers of each participant's sites, and the pairs it fires itself
        marks = {combo: tuple(frozenset(_marker(k) for k in kinds[ci][p][combo[p]]) for p in (0, 1))
                 for combo in combos} if parts == 2 else {}
        blocks = {combo: _unshared_blocks(*m) for combo, m in marks.items()}
        if any(len(b) > 1 for b in blocks.values()):
            kmc_only.append({"template": t.name, "family": t.family, "reason": "deferred_not_a_product"})
            continue
        placed = deferred = False
        for combo in combos:
            markers = None
            if parts == 2:
                if not blocks[combo]:
                    deferred = True        # every pair is a transition of a template it defers to
                    continue
                (markers,) = blocks[combo]
                if markers == marks[combo]:
                    markers = None         # no pair is
                pairs = [(k0, k1) for k0 in kinds[ci][0][combo[0]] for k1 in kinds[ci][1][combo[1]]
                         if markers is None or (_marker(k0) in markers[0] and _marker(k1) in markers[1])]
                leaking = [me._pair_leak(cls, k0, k1) for k0, k1 in pairs]
                if all(leaking):
                    continue
                partial = any(leaking)
            else:
                partial = False
            role = roles(t, ktype)
            order = sorted(range(parts), key=lambda p: (role[p], combo[p]))
            kid = f"{ktype}:" + "+".join(f"{role[p]}={combo[p]}" for p in order)
            k = kernels.setdefault(kid, {"type": ktype, "roles": [role[p] for p in order],
                                         "classes": [combo[p] for p in order],
                                         "molecularity": parts, "templates": []})
            term = {"template": t.name, "phi_thinned": row["phi_thinned"]}
            if parts == 1:
                c = row["participants"][0][combo[0]]
                if combo[0] in chains:
                    term.update({"k_mu1": t.weight * c["a"], "k_mu0": t.weight * c["b"]})
                else:
                    term["k_per_feature"] = t.weight * c["per_feature"]
            else:
                sites = [row["participants"][p][combo[p]] if markers is None
                         else chain_count(ci, p, markers[p], combo[p]) if combo[p] in chains
                         else feature_count(combo[p], ci, p, markers[p]) for p in (0, 1)]
                term.update({"w_over_divisor": t.weight / cls.divisor,
                             "S0": sites[order[0]], "S1": sites[order[1]]})
                if partial:
                    term["partially_leaking"] = True
                if markers is not None:
                    term["sites_deferred_to"] = sorted({templates[c].name for p in (0, 1)
                                                        for m in marks[combo][p] for c in m})
            k["templates"].append(term)
            jobs.append((term, ci, combo, order, markers))
            row["kernels"].append(kid)
            placed = True
        if not placed:
            kmc_only.append({"template": t.name, "family": t.family,
                             "reason": "deferred" if deferred else "branching"})
    for k in kernels.values():
        if k["molecularity"] == 1 and k["classes"][0] in chains:
            k["k_mu1"] = math.fsum(x["k_mu1"] for x in k["templates"])
            k["k_mu0"] = math.fsum(x["k_mu0"] for x in k["templates"])
            k["rate_law"] = "R = k_mu1*mu1 + k_mu0*mu0 (mol/s; mu in mol, n in primitive periods)"
        elif k["molecularity"] == 1:
            k["k_per_feature"] = math.fsum(x["k_per_feature"] for x in k["templates"])
            k["rate_law"] = "R = k_per_feature*mu0(pool carrying the feature) (mol/s)"
        else:
            k["rate_law"] = "R = sum (w_over_divisor)*S0*S1/V (mol/s; S_p in mol of sites, V in m3)"
    unit_mass = spec.repeat_unit.composition.molar_mass_g_mol() / ppu
    laws = _Laws(engine, probes, n_hi)
    measured, placed_overrides = {}, set()
    for i, (term, ci, combo, order, markers) in enumerate(jobs, 1):
        t = templates[ci]
        start = time.perf_counter()
        _law_term(laws, engine, probes, t, term, ci, combo, order, measured, placed_overrides, entries,
                  template_rows[ci]["participants"], markers)
        status = "law" if term["law"] is not None else "inert" if term.get("law_inert") else "refused"
        _progress(f"law {i}/{len(jobs)} {t.name} {'+'.join(combo)} {status} "
                  f"{time.perf_counter() - start:.1f}s")
    unplaced = [dict(o, template=templates[ci].name) for ci, rows in sorted(measured.items()) for o in rows
                if (templates[ci].name, o["name"]) not in placed_overrides]
    families = {}
    for family in FAMILIES:
        labels = {}
        for n in range(1, n_hi + 1):
            try:
                labels[n] = probes.label(family, n)
            except _OffDeclared:
                continue
        offsets = probes.offsets.get(family, set())
        if len(offsets) > 1:
            raise KernelRatesError(f"{family}: the molar mass is not affine in n ({sorted(offsets)})")
        compositions = probes.composition_offsets.get(family, set())
        if len(compositions) > 1:
            raise KernelRatesError(f"{family}: the composition is not affine in n ({sorted(compositions)})")
        if not labels:
            continue
        far = labels[max(labels)]
        families[family] = {
            "label_far": far, "min_n": min(labels),
            "labels": {str(n): v for n, v in sorted(labels.items()) if v != far},
            "molar_mass_offset_g_mol": float(next(iter(offsets))),
            "composition_offset": {el: str(v) for el, v in next(iter(compositions))},
            "declared_species": {f"{lab}@{n}": name for (lab, n), name in sorted(probes.species.items())
                                 if family_of(lab) == family and type_of(lab) is None}}
    extra = {} if probes.legacy else {"chain_types": _chain_types_doc(probes, n_hi),
                                      "unreached": _unreached(kernels, probes)}
    out = {"schema": SCHEMA, "temperature_K": float(temperature),
           "catalog": {"chemistry_sha256": (doc.get("provenance") or {}).get("chemistry_sha256"),
                       "templates": len(templates)},
           "moment_unit": {"periods_per_declared_unit": ppu, "molar_mass_g_mol": unit_mass,
                           "backbone_atoms_per_period": nbb,
                           "composition": {el: str(v) for el, v in sorted(probes.composition.items()) if v},
                           "definition": "n counts primitive periods of the declared backbone"},
           "chain_lengths_probed": list(range(1, n_hi + 1)),
           "families": families, "laws_exact_to_n": n_hi, "overrides_unplaced": unplaced,
           "kernels": dict(sorted(kernels.items())), "kmc_only": kmc_only, "templates": template_rows, **extra}
    if run is not None:
        out["hazard"] = hazard(doc, out, run)
    return out


def template_rows_sha256(table: dict, names: Sequence[str]) -> str:
    """The digest of what a table states for the templates ``names``: each
    one's row (weight, divisor, participants) and its kernel term (product
    law, overrides), canonical JSON in name order -- the same in a copy of
    the table cut to fewer kernels or templates."""
    rows = {r["template"]: r for r in table["templates"]}
    terms = {t["template"]: dict(t, kernel=kid) for kid, k in table["kernels"].items() for t in k["templates"]}
    missing = sorted(set(names) - (set(rows) & set(terms)))
    if missing:
        raise KernelRatesError(f"templates {missing} are not in the table")
    body = [[rows[n], terms[n]] for n in sorted(set(names))]
    return hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def hazard(doc: dict, table: dict, run: dict) -> dict:
    """Fired counts (gross and netted over reciprocal pairs) and the pooled
    leak fraction of a ``kmc-melt-run/1`` manifest, by template, kernel and
    kMC-only entry."""
    fired = Counter()
    for r in run.get("runs", []):
        fired.update((r.get("report") or {}).get("template_counts", {}))
    total = sum(fired.values())
    by_kernel = Counter()
    for row in table["templates"]:
        if len(row["kernels"]) == 1:
            by_kernel[row["kernels"][0]] += fired.get(row["template"], 0)
        elif row["kernels"]:
            by_kernel["several:" + "|".join(row["kernels"])] += fired.get(row["template"], 0)
    only = [{"template": e["template"], "reason": e["reason"], "fired": fired.get(e["template"], 0),
             "share": fired.get(e["template"], 0) / total if total else None} for e in table["kmc_only"]]
    record = Counter()
    for name, v in fired.items():
        record[name.split(":")[0]] += v
    pairs, seen = [], set()
    for pay in doc["payloads"]:
        rec = pay["record_event_id"][4:16]
        rev = (pay.get("reverse_of") or "")[4:16]
        if not rev or rec in seen:
            continue
        seen.update({rec, rev})
        pairs.append({"records": [rec, rev], "fired": [record[rec], record[rev]],
                      "net": abs(record[rec] - record[rev])})
    leak = (run.get("leak") or {}).get("pooled") or {}
    return {"events": total, "by_kernel_gross": dict(by_kernel), "kmc_only": only,
            "kmc_only_share_gross": math.fsum(x["fired"] for x in only) / total if total else None,
            "pairs_net": pairs, "net_events": sum(p["net"] for p in pairs),
            "leak_fraction": leak.get("fraction"),
            "note": "leak hazard never fires; it is outside every kernel and every count above"}


def main(argv=None) -> int:
    from rmgpy.kmc import run as kr
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.kernel_rates",
                                     description="MoM kernel rates from a kMC template catalog")
    parser.add_argument("deck")
    parser.add_argument("catalog")
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--run", default=None, help="kmc-melt-run/1 run.json for hazard shares")
    parser.add_argument("--out", default=None, help="default: <catalog>.kernel_rates.json beside it")
    parser.add_argument("--allow-unbound-catalog", action="store_true")
    args = parser.parse_args(argv)
    try:
        deck = kr.read_deck(args.deck)
        doc = json.loads(Path(args.catalog).read_text())
        run = json.loads(Path(args.run).read_text()) if args.run else None
        table = kernel_rates(doc, deck.spec, args.temperature, run=run,
                             allow_unbound=args.allow_unbound_catalog)
    except (KernelRatesError, kr.RunError, OSError, ValueError) as error:
        print(f"refused: {error}", file=sys.stderr)
        return 2
    out = Path(args.out) if args.out else Path(args.catalog).with_suffix(".kernel_rates.json")
    out.write_text(json.dumps(table, indent=1, sort_keys=True) + "\n")
    print(f"written: {out} ({len(table['kernels'])} kernels, {len(table['kmc_only'])} kMC-only)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
