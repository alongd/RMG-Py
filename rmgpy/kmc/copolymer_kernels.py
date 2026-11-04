"""MoM kernel rates of a linear statistical copolymer from a kMC template
catalog (I-240, slice S5 of the copolymer design).

``copolymer_kernel_rates(doc, spec, temperature)`` places the templates a
``kmc-template-catalog/1`` loads at ``temperature`` on one copolymer pool
(``rmgpy.copolymer_pool``: one residue per unit, moments counting residues)
or refuses them by name.  Nothing is fitted and nothing is dropped.

**Period.**  One backbone atom (``backbone_atoms_per_period`` 1); every
residue holds the same number of backbone atoms (the seed refuses otherwise),
which is ``periods_per_declared_unit``.  ``N`` counts periods.

**Weighting rule.**  A template's sites on a chain of n residues are its
expected count over every residue word of length n, each at its frequency in
the declared sequence statistics (``declared_environments.word_frequencies``:
pi(w1) T(w1, w2) ... T(w_{n-1}, w_n), the chance that a chain of n residues
is that word), each word's count read by the melt engine's own index
(``MeltEngine._scan``) on that word materialized as one chain
(``melt_state.materialize_sequences``)::

    s_T(n) = sum_w f(w) c_T(w)

For every n probed (1 .. ``units_hi``) this is exact in rationals; the row's
``sites`` are ``a N + b`` (N = n * periods per residue), required to hold
exactly on the trailing ``EXACT_POINTS`` lengths (refused otherwise) and from
``exact_from_n`` periods on.  Each window of the chain is thereby counted at
its own frequency: a site in a window of residues u is weighted by f(u).

**Kernel.**  Chain homolysis (``kernel_rates.kernel_type``) of a backbone
bond is the one kernel, ``homolysis:substrate=chain``.  Its MoM law is the
uniform cut: the expected propensity at each cut position j (fragments of j
and N - j periods), ``r_N(j) = sum_T w_T sum_w f(w) c_T(w, j)``, must be one
value A on min_fragment <= j <= N - min_fragment and 0 elsewhere, on every
probed length from the shortest initial chain (or the trailing lengths when
every initial chain is longer than the probe), and A must be the same at
every length (refused otherwise, ``refused:nonuniform_cut``).  Templates are
grouped by ``phi_thinned`` (each group is checked on its own) and each group
is one ``kernel_homolysis`` term (``mom_homolysis_terms``): one site per cut
position (a = 1, b = 1 - 2 min_fragment) at rate w = A.

**Drift.**  ``per_residue_rate`` (1/s per residue) is a residue's marginal
propensity: ``ppu * a + E[R | s_i = r] - E[R]``, with R = sum_T w_T c_T the
chain's homolysis propensity and i the middle residue of the longest probed
chain (exact for statistics without correlation beyond the windows, e.g.
``ideal_random``; for ``markov1`` the value at that length of a geometric
series).  The kernel is ``residue_selective`` when these differ; the
selective part of the propensity is what lies above the slowest residue's
rate, ``a_floor = min_r per_residue_rate[r] / ppu`` per period, and its
fraction ``(a - a_floor) mu1 / (a mu1 + b mu0)`` (mu in periods;
``composition_drift``) is reported under ``limited:composition_drift``:
one pool decays every chain at the mean rate, while the chains that survive
selective chemistry are richer in the slower residue and decay slower (the
sequence-level decay), and the composition a residue's loss from the polymer
moves is a kMC observable that MoM withholds.  Never corrected.

**Scope.**  Like kernel mode's homolysis on a homopolymer, the kernel cuts
the chain pool only: a fragment in the end-radical pool is not cut again
(its cut would make a two-radical fragment, no pool of this mode), so MoM's
chain count stays below twice its initial value while kMC's grows past it.

**Refused** (``refused``, with the reason): a template with no kernel type
(``no_kernel``); one with two participants (``refused:residue_context``: its
sites are features -- radical ends, interior radicals, alkenes -- or pairs
of them, each on one residue or another, which one pool does not carry); one
with no site on a pristine chain of any declared word
(``no_site_on_declared_states``, the residue context named); a chain template whose transitions
leak (``branching``); a chain template of another type, or a homolysis of a
bond off the backbone (``refused:no_copolymer_kernel``).  The kernel is
refused whole (``refused:nonuniform_cut``) when its cut law is not uniform.
"""

import math
from collections import Counter, defaultdict
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

import numpy

from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_state as ms
from rmgpy.kmc.declared_environments import word_frequencies

SCHEMA = "kmc-copolymer-kernel-rates/1"
EXACT_POINTS = K.EXACT_POINTS
COMPOSITION_DRIFT = "limited:composition_drift"
RESIDUE_CONTEXT = "refused:residue_context"
NO_COPOLYMER_KERNEL = "refused:no_copolymer_kernel"
NONUNIFORM_CUT = "refused:nonuniform_cut"
HOMOLYSIS = "homolysis:substrate=chain"
MAX_CHAINS = 5000          # residue words materialized at once (sum over the probed lengths)


class CopolymerKernelsError(ValueError):
    """The catalog, the declaration or a count is refused."""


def _backbone_path(state, cid) -> List[int]:
    """The chain's backbone atoms from its head end (unit 1)."""
    origin = state.origin
    bb = {a for a in state.components.atoms(cid) if origin[a].kind == "chain" and origin[a].is_backbone}
    nbrs = {a: [b for b in state.atoms[a].edges if b in bb] for a in bb}
    ends = [a for a in bb if len(nbrs[a]) <= 1]
    head = min(ends, key=lambda a: int(origin[a].segment.split()[1]))
    path, prev = [head], None
    while len(path) < len(bb):
        nxt = [b for b in nbrs[path[-1]] if b != prev]
        if len(nxt) != 1:
            raise CopolymerKernelsError("a declared chain's backbone is not a simple path")
        prev = path[-1]
        path.append(nxt[0])
    return path


def _affine(series: Dict[int, Fraction]) -> dict:
    """a, b of s(N) = a N + b exact on the trailing EXACT_POINTS lengths, and
    the smallest N from which it holds; refused when not exact."""
    pts = sorted(series.items())
    if len(pts) < EXACT_POINTS:
        raise CopolymerKernelsError("too few chain lengths to establish a site polynomial")
    (n1, y1), (n2, y2) = pts[-2], pts[-1]
    a = Fraction(y2 - y1, n2 - n1)
    b = y1 - a * n1
    if any(y != a * n + b for n, y in pts[-EXACT_POINTS:]):
        raise CopolymerKernelsError(f"expected site counts {dict(pts)} are not affine in N on the trailing lengths")
    start = pts[-EXACT_POINTS][0]
    for n, y in reversed(pts[:-EXACT_POINTS]):
        if y != a * n + b:
            break
        start = n
    return {"a": a, "b": b, "exact_from_n": start}


def _uniform(profiles: Dict[int, Dict[int, Fraction]], lengths: Sequence[int], where: str) -> Tuple[Fraction, int]:
    """(A, e) of profiles r_N(j) = A on e <= j <= N - e, 0 elsewhere, one A and
    one e on every given length N; refused otherwise."""
    found = set()
    for N in lengths:
        profile = profiles[N]
        live = sorted(j for j, v in profile.items() if v)
        if not live:
            raise CopolymerKernelsError(f"{NONUNIFORM_CUT}: {where} cuts no chain of {N} periods")
        e = live[0]
        values = {profile.get(j, Fraction(0)) for j in range(e, N - e + 1)}
        if live != list(range(e, N - e + 1)) or len(values) != 1:
            raise CopolymerKernelsError(
                f"{NONUNIFORM_CUT}: {where}: the expected propensity per cut position on chains of {N} periods "
                f"{dict(sorted(profile.items()))} is not one value on min_fragment..N - min_fragment")
        found.add((values.pop(), e))
    if len(found) != 1:
        raise CopolymerKernelsError(f"{NONUNIFORM_CUT}: {where}: the uniform cut differs between chain lengths "
                                    f"{sorted(found)}")
    return found.pop()


def copolymer_kernel_rates(doc: dict, spec, temperature: float, *, allow_unbound: bool = False,
                           max_units: Optional[int] = None) -> dict:
    """The copolymer kernel-rate table (module doc) of a catalog document
    for one ``PolymerSpec`` at one temperature."""
    from rmgpy.copolymer_pool import copolymer_pool_seed
    from rmgpy.kmc import run as kr
    from rmgpy.polymer_spec import PolymerSpec
    if not isinstance(spec, PolymerSpec):
        raise CopolymerKernelsError("copolymer kernel rates need a PolymerSpec (schema v2); a linear homopolymer "
                                    "declaration is kernel_rates.kernel_rates'")
    if doc.get("schema") != "kmc-template-catalog/1" or not isinstance(doc.get("payloads"), list):
        raise CopolymerKernelsError("not a kmc-template-catalog/1 document")
    bound = (doc.get("provenance") or {}).get("chemistry_sha256")
    if bound != kr.chemistry_sha256(spec) and not allow_unbound:
        raise CopolymerKernelsError(f"the catalog's chemistry digest {bound!r} is not the declaration's "
                                    f"({kr.chemistry_sha256(spec)})")
    seed = copolymer_pool_seed(spec)
    ppu = seed.backbone_atoms_per_residue
    catalog = me.PayloadCatalog(doc["payloads"])
    templates = list(catalog.templates(float(temperature)))
    window = max(sum(1 for a in t.atoms if a.element != "H") for t in templates) if templates else 1
    units_hi = max_units or 2 * math.ceil(window / ppu) + EXACT_POINTS
    words = {n: sorted(word_frequencies(spec, n).items()) for n in range(1, units_hi + 1)}
    if sum(len(w) for w in words.values()) > MAX_CHAINS:
        raise CopolymerKernelsError(f"{sum(len(w) for w in words.values())} residue words up to {units_hi} "
                                    f"residues exceed {MAX_CHAINS} chains (max_units)")
    order = [(n, word, f) for n in sorted(words) for word, f in words[n]]
    state = ms.materialize_sequences(spec, [word for _n, word, _f in order])
    engine = me.MeltEngine(spec, catalog, float(temperature), numpy.random.default_rng(0), state=state)
    if [t.name for t in engine.templates] != [t.name for t in templates]:
        raise CopolymerKernelsError("the engine loaded another template set")
    classes = engine.classes
    weight = [Fraction(float(t.weight)) for t in templates]

    # per word: per template its count, its cut positions (homolysis of a backbone bond), its leaks
    mean = defaultdict(lambda: defaultdict(Fraction))         # template -> N -> sum_w f c_T(w)
    cuts = defaultdict(lambda: defaultdict(lambda: defaultdict(Fraction)))   # template -> N -> j -> sum_w f c
    off_backbone, leak, paired = Counter(), Counter(), Counter()
    per_word = []                                             # (n, word, f, {template: count})
    for (n, word, f), cid in zip(order, sorted(state.components.ids())):
        path = _backbone_path(state, cid)
        if len(path) != n * ppu:
            raise CopolymerKernelsError(f"word {word}: {len(path)} backbone atoms, not {n * ppu}")
        position = {a: i + 1 for i, a in enumerate(path)}
        entry = engine._scan([cid], state=state)[cid]
        counts = Counter()
        for ci, (allowed, leaking) in K._uni_all(entry).items():
            leak[ci] += leaking
            for tr in allowed:
                counts[ci] += 1
                broken = [e for e in tr.edits if isinstance(e, ms.SetBond) and e.new == 0
                          and e.a in position and e.b in position and abs(position[e.a] - position[e.b]) == 1]
                if len(broken) != 1:
                    off_backbone[ci] += 1
                    continue
                cuts[ci][n * ppu][min(position[broken[0].a], position[broken[0].b])] += f
        for (ci, _p, _kind), sites in entry.sites.items():
            paired[ci] += len(sites)
        for ci, c in counts.items():
            mean[ci][n * ppu] += f * c
        per_word.append((n, word, f, counts))

    rows, refused, placed = [], [], []
    for ci, (t, cls) in enumerate(zip(templates, classes)):
        ktype = K.kernel_type(t, cls)
        row = {"template": t.name, "family": t.family, "type": ktype, "weight": float(t.weight),
               "weight_units": "m3/(mol*s) per matched transition" if cls.bimolecular else
               "1/s per matched transition", "phi_thinned": bool(cls.met_reverse), "kernels": []}
        series = {n * ppu: mean[ci][n * ppu] for n in words}
        fit = _affine(series)
        row["sites"] = {"a": str(fit["a"]), "b": str(fit["b"]), "exact_from_n": fit["exact_from_n"]}
        rows.append(row)
        reason = None
        if ktype is None:
            reason = "no_kernel"
        elif cls.bimolecular or len(t.participants) != 1:
            reason = (f"{RESIDUE_CONTEXT}: a two-participant {ktype} pairs feature sites (radical ends, interior "
                      "radicals, alkenes), each on one residue or another; one copolymer pool carries no "
                      "residue context")
        elif not any(series.values()):
            reason = ("no_site_on_declared_states: no site on a pristine chain of any declared word; a site on a "
                      "feature (radical end, interior radical, alkene) would sit on one residue or another, "
                      f"which one copolymer pool does not carry ({RESIDUE_CONTEXT})")
        elif leak[ci]:
            reason = "branching"
        elif ktype != "homolysis" or off_backbone[ci]:
            reason = (f"{NO_COPOLYMER_KERNEL}: only the homolysis of a backbone bond has a one-pool form "
                      f"(type {ktype}, {off_backbone[ci]} cuts off the backbone)")
        if reason is not None:
            refused.append({"template": t.name, "family": t.family, "reason": reason})
        else:
            placed.append(ci)

    kernels, limitations = {}, []
    if placed:
        n_min = min(dp for dp, _c in seed.dp_counts) * ppu
        lengths = sorted(n * ppu for n in words)
        checked = [N for N in lengths if N >= n_min] or lengths[-EXACT_POINTS:]
        groups = []
        for thinned in (False, True):
            members = [ci for ci in placed if bool(classes[ci].met_reverse) == thinned]
            if not members:
                continue
            profiles = {N: defaultdict(Fraction) for N in lengths}
            for ci in members:
                for N, by_j in cuts[ci].items():
                    for j, v in by_j.items():
                        profiles[N][j] += weight[ci] * v
            where = f"{HOMOLYSIS} ({'phi_thinned' if thinned else 'unthinned'} templates)"
            try:
                A, e = _uniform(profiles, checked, where)
            except CopolymerKernelsError as error:
                refused.extend({"template": templates[ci].name, "family": templates[ci].family,
                                "reason": str(error)} for ci in members)
                continue
            if n_min < 2 * e - 1:
                raise CopolymerKernelsError(f"{where}: initial chains of {n_min} periods are shorter than "
                                            f"{2 * e - 1}, where the cut law counts negative sites")
            groups.append({"phi_thinned": thinned, "rate_per_position": A, "min_fragment": e,
                           "templates": [templates[ci].name for ci in members], "members": members})
        if groups:
            if len({g["min_fragment"] for g in groups}) != 1:
                raise CopolymerKernelsError(f"{HOMOLYSIS}: thinned and unthinned templates cut with different "
                                            f"min_fragment {sorted({g['min_fragment'] for g in groups})}")
            members = [ci for g in groups for ci in g["members"]]
            e = groups[0]["min_fragment"]
            a = sum((g["rate_per_position"] for g in groups), Fraction(0))
            b = a * (1 - 2 * e)
            # the marginal propensity of a residue at the middle of the longest probed chain
            n_hi = max(words)
            middle = (n_hi - 1) // 2
            total, given, mass = Fraction(0), defaultdict(Fraction), defaultdict(Fraction)
            for n, word, f, counts in per_word:
                if n != n_hi:
                    continue
                R = sum((weight[ci] * counts[ci] for ci in members), Fraction(0))
                total += f * R
                given[word[middle]] += f * R
                mass[word[middle]] += f
            rates = {r: ppu * a + given[r] / mass[r] - total for r in sorted(mass)}
            floor = min(rates.values()) / ppu
            selective = len(set(rates.values())) > 1
            for ci in members:
                rows[ci]["kernels"].append(HOMOLYSIS)
            names = [templates[ci].name for ci in members]
            kernels[HOMOLYSIS] = {
                "type": "homolysis", "roles": ["substrate"], "classes": ["chain"], "molecularity": 1,
                "templates": [{"template": rows[ci]["template"], "phi_thinned": rows[ci]["phi_thinned"]}
                              for ci in members],
                "rate": {"a": float(a), "b": float(b)},
                "rate_law": "R = a*N + b per chain (1/s; N in periods, intrinsic: phi not applied)",
                "mom": {"min_fragment": e, "valid_from_periods": 2 * e - 1,
                        "terms": [{"phi_thinned": g["phi_thinned"], "w": str(g["rate_per_position"]),
                                   "templates": g["templates"]} for g in groups]},
                "residue_selective": selective,
                "per_residue_rate": {r: float(v) for r, v in rates.items()},
                "drift": {"a": str(a), "b": str(b), "a_floor": str(floor)}}
            if selective:
                limitations.append({
                    "name": COMPOSITION_DRIFT, "kernels": [HOMOLYSIS], "templates": names,
                    "per_residue_rate": {r: float(v) for r, v in rates.items()},
                    "definition": "selective fraction = (a - a_floor) mu1 / (a mu1 + b mu0), mu in periods, "
                                  "a_floor = min_r per_residue_rate[r] / periods_per_declared_unit; the "
                                  "composition this drifts is withheld from MoM, never corrected"})
    return {"schema": SCHEMA, "temperature_K": float(temperature),
            "catalog": {"chemistry_sha256": bound, "templates": len(templates)},
            "moment_unit": {"periods_per_declared_unit": ppu, "backbone_atoms_per_period": 1,
                            "molar_mass_g_mol": float(seed.mean_residue_molar_mass) / ppu,
                            "definition": "n counts backbone atoms; a declared unit is one residue of mean "
                                          "molar mass sum x_i M_i"},
            "weighting": "s_T(n) = sum over residue words w of length n of word_frequencies(w) * c_T(w)",
            "residue_words_probed": {str(n): len(w) for n, w in words.items()},
            "kernels": kernels, "refused": refused, "limitations": limitations, "templates": rows}


def mom_homolysis_terms(table: dict, *, phi: float, daughter: str) -> List[dict]:
    """polymer-mom's ``kernel_homolysis`` terms of the table's chain
    homolysis: per phi_thinned group one site per cut position (a = 1,
    b = 1 - 2 min_fragment) at the group's rate per position w, phi on the
    thinned group."""
    if table.get("schema") != SCHEMA:
        raise CopolymerKernelsError(f"not a {SCHEMA} table (got {table.get('schema')!r})")
    kernel = table["kernels"].get(HOMOLYSIS)
    if kernel is None:
        raise CopolymerKernelsError(f"the table has no {HOMOLYSIS} kernel")
    e = int(kernel["mom"]["min_fragment"])
    ppu = int(table["moment_unit"]["periods_per_declared_unit"])
    terms = kernel["mom"]["terms"]
    return [{"template": HOMOLYSIS if len(terms) == 1 else
             f"{HOMOLYSIS}|{'phi_thinned' if g['phi_thinned'] else 'unthinned'}",
             "w": float(Fraction(g["w"])), "phi": float(phi) if g["phi_thinned"] else 1.0,
             "a": 1.0, "b": float(1 - 2 * e), "min_fragment": e, "periods_per_unit": ppu,
             "daughter": daughter, "valid_from_periods": 2 * e - 1} for g in terms]


def composition_drift(table: dict, samples: Sequence[Tuple[float, Tuple[float, float]]]) -> dict:
    """``limited:composition_drift`` over (time, (mu0, mu1)) of the chain pool
    (mu in residues): the selective fraction of the homolysis propensity at
    each time (module doc)."""
    ppu = int(table["moment_unit"]["periods_per_declared_unit"])
    series = []
    for time, (mu0, mu1) in samples:
        above = total = 0.0
        for kernel in table["kernels"].values():
            d = kernel["drift"]
            a, b, floor = float(Fraction(d["a"])), float(Fraction(d["b"])), float(Fraction(d["a_floor"]))
            above += (a - floor) * ppu * mu1
            total += a * ppu * mu1 + b * mu0
        series.append({"time_s": float(time), "selective_fraction": above / total if total > 0 else None})
    return {"name": COMPOSITION_DRIFT, "series": series}
