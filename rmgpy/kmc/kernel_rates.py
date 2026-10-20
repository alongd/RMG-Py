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
  branching kind) on chains of every length up to ``n_max`` and required to be
  exactly affine on the trailing lengths (refused otherwise: never fitted);
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
with no site on any state built here, e.g. a homolysis at a branch point).

With a ``kmc-melt-run/1`` manifest, ``hazard`` gives the fired counts per
template, kernel and kMC-only entry (gross), the counts netted over each
reciprocal record pair (``reverse_of``), and the run's pooled leak fraction.
"""

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy

from rmgpy.kmc import melt_engine as me
from rmgpy.kmc.declared_chain import ChainGeometry
from rmgpy.kmc.melt_state import SetBond, SetRadical, materialize

SCHEMA = "kmc-kernel-rates/1"
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
        for ci, (allowed, leak) in entry.uni.items():
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
                out[(ci, p, kind)] += n
        return out


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


def _states(spec, n_values, counter, ppu: int):
    """``chain``: {(class index, participant, kind): {n: count}} on pristine
    chains (n in primitive periods); ``features``: {class label: {"far":
    {units: counts}, "profile": {position: [counts, ...]}}} (``_features``)."""
    chain: Dict[Tuple, Dict[int, int]] = defaultdict(dict)
    features: Dict[str, dict] = {}
    longest = max(n_values)
    for units in n_values:
        n = units * ppu
        state = materialize(_single_chain_spec(spec, units))
        path = _backbone_order(state)
        cid = state.components.component_of(path[0])
        for key, v in counter.counts(state, cid).items():
            chain[key][n] = v
        for key in chain:
            chain[key].setdefault(n, 0)
        _features(spec, units, counter, features, far_ok=(units > longest - EXACT_POINTS))
    return chain, features


def _features(spec, units, counter, features, far_ok):
    """Feature counts on one chain of ``units`` declared units, each feature
    on a fresh state, at every position: ``("fragment_backbone_atoms", m)``
    for an end feature on a fragment of m backbone atoms, and
    ``("chain_backbone_atoms", nb, "end_distance", d)`` for an interior feature
    d backbone bonds from the nearer end of a chain of nb backbone atoms.
    With ``far_ok`` the features in the middle of the chain (a cut of a
    middle bond; the middle backbone atoms or bonds) also give the far-field
    counts."""
    probe = materialize(_single_chain_spec(spec, units))
    nb = len(_backbone_order(probe))
    per_unit = max(1, nb // units)
    middle = range((units // 2) * per_unit, (units // 2) * per_unit + per_unit)

    def fresh():
        s = materialize(_single_chain_spec(spec, units))
        return s, _backbone_order(s)

    def note(kind, label, counts, position, far):
        slot = features.setdefault(f"{kind}:{label}", {"far": {}, "profile": {}})
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
                r, q, m = u, (p[k - 2] if k >= 2 else None), k
            else:
                r, q, m = v, (p[k + 1] if k + 1 < nb else None), nb - k
            far = (k - 1) in middle
            position = ("fragment_backbone_atoms", m)
            note("radical_end", _label(s, [r]), counter.counts(s, s.components.component_of(r), frozenset([r])),
                 position, far)
            if q is not None and _hydrogens(s, q):
                hq = _hydrogens(s, q)[0]
                s.apply([SetBond(q, hq, 1.0, 0.0), SetRadical(hq, 0, 1), SetBond(r, q, 1.0, 2.0),
                         SetRadical(r, s.atoms[r].radical, s.atoms[r].radical - 1)])
                note("alkene_terminal", _label(s, [r, q]),
                     counter.counts(s, s.components.component_of(r), frozenset([r, q])), position, far)

    for j in range(1, nb - 1):
        far = j in middle
        position = ("chain_backbone_atoms", nb, "end_distance", min(j, nb - 1 - j))
        s, p = fresh()
        edits = _remove_h(s, p[j])
        if edits:
            s.apply(edits)
            note("radical_interior", _label(s, [p[j]]),
                 counter.counts(s, s.components.component_of(p[j]), frozenset([p[j]])), position, far)
        if j + 1 < nb - 1:
            s, p = fresh()
            a, b = p[j], p[j + 1]
            ha, hb = _hydrogens(s, a), _hydrogens(s, b)
            if ha and hb and s.atoms[a].edges.get(b) == 1.0:
                s.apply([SetBond(a, ha[0], 1.0, 0.0), SetBond(b, hb[0], 1.0, 0.0),
                         SetRadical(ha[0], 0, 1), SetRadical(hb[0], 0, 1), SetBond(a, b, 1.0, 2.0)])
                note("alkene_internal", _label(s, [a, b]),
                     counter.counts(s, s.components.component_of(a), frozenset([a, b])), position, far)


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
    window = max(sum(1 for a in t.atoms if a.element != "H") for t in templates) if templates else 1
    bb_per_unit = len(geometry.backbone_keys)
    units_hi = max_units or max(2 * EXACT_POINTS + 2, math.ceil((2 * window + EXTRA_UNITS) / bb_per_unit) + 1)
    engine = me.MeltEngine(_single_chain_spec(spec, units_hi), catalog, float(temperature),
                           numpy.random.default_rng(0))
    if [t.name for t in engine.templates] != [t.name for t in templates]:
        raise KernelRatesError("the engine loaded another template set")
    counter = _Counter(engine)
    n_values = list(range(1, units_hi + 1))
    chain, features = _states(spec, n_values, counter, ppu)
    classes = engine.classes

    # per template, participant -> {class: count description}
    by_template = {ci: defaultdict(dict) for ci in range(len(templates))}
    kinds = {ci: defaultdict(lambda: defaultdict(set)) for ci in range(len(templates))}
    for (ci, p, kind), series in chain.items():
        if any(series.values()):
            by_template[ci][p].setdefault("chain", Counter())
            kinds[ci][p]["chain"].add(kind)
    chain_poly = {}
    for ci in range(len(templates)):
        for p in list(by_template[ci]):
            if "chain" in by_template[ci][p]:
                total = Counter()
                for (cj, q, kind), series in chain.items():
                    if cj == ci and q == p:
                        for n, v in series.items():
                            total[n] += v
                chain_poly[(ci, p)] = _affine(dict(total))
    feature_counts = {}
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
        keys = set(far)
        for counts in slot["profile"].values():
            for c in counts:
                keys.update(c)
        for ci, p in sorted({(k[0], k[1]) for k in keys}):
            for k in keys:
                if k[0] == ci and k[1] == p:
                    kinds[ci][p][name].add(k[2])
            c_far = sum(v for k, v in far.items() if k[0] == ci and k[1] == p)
            prof = {}
            for position, counts in sorted(slot["profile"].items()):
                values = sorted({sum(v for k, v in c.items() if k[0] == ci and k[1] == p) for c in counts})
                if values != [c_far]:
                    prof[",".join(f"{position[i]}={position[i + 1]}" for i in range(0, len(position), 2))] = values
            if c_far or any(any(v) for v in prof.values()):
                feature_counts[(ci, p, name)] = {"per_feature": c_far, "profile": prof}
                by_template[ci][p][name] = Counter()

    kernels: Dict[str, dict] = {}
    kmc_only: List[dict] = []
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
                if name == "chain":
                    entry[name] = chain_poly[(ci, p)]
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
        placed = False
        for combo in combos:
            if parts == 2:
                pairs = [(k0, k1) for k0 in kinds[ci][0][combo[0]] for k1 in kinds[ci][1][combo[1]]]
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
                if combo[0] == "chain":
                    term.update({"k_mu1": t.weight * c["a"], "k_mu0": t.weight * c["b"]})
                else:
                    term["k_per_feature"] = t.weight * c["per_feature"]
            else:
                term.update({"w_over_divisor": t.weight / cls.divisor,
                             "S0": row["participants"][order[0]][combo[order[0]]],
                             "S1": row["participants"][order[1]][combo[order[1]]]})
                if partial:
                    term["partially_leaking"] = True
            k["templates"].append(term)
            row["kernels"].append(kid)
            placed = True
        if not placed:
            kmc_only.append({"template": t.name, "family": t.family, "reason": "branching"})
    for k in kernels.values():
        if k["molecularity"] == 1 and k["classes"] == ["chain"]:
            k["k_mu1"] = math.fsum(x["k_mu1"] for x in k["templates"])
            k["k_mu0"] = math.fsum(x["k_mu0"] for x in k["templates"])
            k["rate_law"] = "R = k_mu1*mu1 + k_mu0*mu0 (mol/s; mu in mol, n in primitive periods)"
        elif k["molecularity"] == 1:
            k["k_per_feature"] = math.fsum(x["k_per_feature"] for x in k["templates"])
            k["rate_law"] = "R = k_per_feature*mu0(pool carrying the feature) (mol/s)"
        else:
            k["rate_law"] = "R = sum (w_over_divisor)*S0*S1/V (mol/s; S_p in mol of sites, V in m3)"
    unit_mass = spec.repeat_unit.composition.molar_mass_g_mol() / ppu
    out = {"schema": SCHEMA, "temperature_K": float(temperature),
           "catalog": {"chemistry_sha256": (doc.get("provenance") or {}).get("chemistry_sha256"),
                       "templates": len(templates)},
           "moment_unit": {"periods_per_declared_unit": ppu, "molar_mass_g_mol": unit_mass,
                           "definition": "n counts primitive periods of the declared backbone"},
           "chain_lengths_probed": [n * ppu for n in n_values],
           "kernels": dict(sorted(kernels.items())), "kmc_only": kmc_only, "templates": template_rows}
    if run is not None:
        out["hazard"] = hazard(doc, out, run)
    return out


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
