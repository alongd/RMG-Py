"""Kernel references (I-214, slice S2 of the shared-chemistry design): an
exact length-resolved population balance (V0) and a length-level stochastic
simulation (V1) of one kernel set, both read from a kernel-rates table
(``kernel_rates``, ``kmc-kernel-rates/3``, or a ``/4`` table in the ``/3``
state model: one chain type, its own reverse; a table of more chain types is
refused by name) and a deck's initial ensemble.

**The kernel model.** A molecule is a state (family, n): its pool family
(``chain``, ``radical_end``, ``alkene_terminal``, ``radical_interior``,
``alkene_internal``) and its length n in primitive periods; its class label
is the table's (``families``: a short fragment may carry its own, e.g. the
methyl radical).  Nothing inside a molecule is resolved.  Each kernel term
of the table acts on the states whose label is its participant class, at
the template weight ``w`` (times ``phi`` on a radical-creating term):

* chain and end families: the term's exact product law at every n the law
  was measured on (``exact``, every site with its products) and its far rule
  beyond (``kernel_rates.predict``);
* interior families: the far-field outcomes (``per_feature`` sites), the
  feature at a position uniform over the positions where an outcome's
  products exist (relative lengths ``alpha*L + beta*x + gamma`` in backbone
  atoms); this position average is the kernel model's, not the melt's;
* unimolecular rate ``w * sites`` per molecule; bimolecular rate
  ``(w / divisor) * S0 * S1 / (N_A * V)`` per volume (S_p the participant's
  sites over the ensemble); a joined product (recombination) is one state of
  length n0 + n1 + gamma;
* a depository override replaces ``w`` by its own weight on its matched
  transitions at its exact sizes (``overrides``).

Terms the table marks inert (an interior participant with no far-field
site: rate 0 by the table's own rate law) and terms whose products are no
declared state (``law_refused``) are listed in every output, never silently
dropped; a kernel filter lists what it excludes.

**PBE** (``population_balance``): the mean-field ODE on state counts of the
deck's representative system (``dN/dt``; bimolecular pairs ``S0*S1``),
backward differentiation with an analytic Jacobian.  Products longer than
``n_max`` periods leave the grid into a per-family overflow (molecules,
periods and periods squared accumulated exactly; they no longer react) and
are reported.

**SSA** (``stochastic``): Gillespie's direct method on the same state
counts; bimolecular propensities ``(w / divisor) * (S0*S1 - X) / (N_A*V)``
with ``X`` the same-molecule pairs (the melt engine's form), pairs drawn by
sites with same-molecule rejection.

**Observables**: the ``whole_*`` keys of ``melt_observables.
mom_reported_quantities`` (``WHOLE_ENSEMBLE_DEFINITION``): mass-DP moments
of every polymer-derived state (not a declared species, at least one
declared unit of backbone), in mol.  For the comparison rows
(``comparison_rows``, I-217) every time also carries ``row``
(``KernelModel.row_observation``): atoms per element (a state holds n *
the period's composition + its family's offset, from the table), mass,
declared species, polymer-derived states per length, overflow molecules.

**Element balance**: at model build every term is checked to conserve each
element exactly (rationals): a unimolecular outcome's products against its
reactant, a join against both reactants, and a separable bimolecular term
by the change of each participant, which must be one change per participant
and sum to none.  A term that does not balance is listed
(``KernelModel.unbalanced``) and both references then show the residual.
"""

import argparse
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy
from scipy import sparse
from scipy.integrate import solve_ivp

from rmgpy import constants
from rmgpy.kmc import kernel_rates as K
from rmgpy.kmc.melt_observables import WHOLE_ENSEMBLE_KEYS, whole_averages

SCHEMA = "kmc-kernel-reference/2"


class KernelReferenceError(ValueError):
    """The table, the deck or a run is refused."""


# --------------------------------------------------------------------------
# the compiled kernel model
# --------------------------------------------------------------------------
class _Bi:
    """One bimolecular term: rate constant (per unordered site pair, over
    N_A*V), per-slot site vectors and per-slot outcomes, or a join."""

    def __init__(self, kernel, template, k, sites, outcomes, join):
        self.kernel, self.template, self.k = kernel, template, k
        self.sites = sites                  # [s0, s1]: numpy arrays over states
        self.outcomes = outcomes            # [ {state: [(mult, target)]} ] per slot (separable)
        self.join = join                    # (family index, gamma in periods) or None


class KernelModel:
    """The kernel model (module doc) of a table for a deck."""

    def __init__(self, table: dict, spec, *, phi: float = 1.0, n_max: Optional[int] = None,
                 kernels: Optional[Sequence[str]] = None, allow_unbound: bool = False):
        from rmgpy.kmc import run as kr
        from rmgpy.kmc.melt_engine import declared_initialization
        from rmgpy.polymer_spec import LinearHomopolymerSpec
        if table.get("schema") != K.SCHEMA:
            raise KernelReferenceError(f"not a {K.SCHEMA} table (got {table.get('schema')!r})")
        beyond = K.state_model_beyond_v3(table)
        if beyond is not None:
            raise KernelReferenceError(f"{beyond}: a state here is (family, n) on that one type")
        if not isinstance(spec, LinearHomopolymerSpec):
            raise KernelReferenceError("the kernel references need a linear homopolymer declaration")
        bound = (table.get("catalog") or {}).get("chemistry_sha256")
        if bound != kr.chemistry_sha256(spec) and not allow_unbound:
            raise KernelReferenceError(f"the table's chemistry digest {bound!r} is not the deck's "
                                       f"({kr.chemistry_sha256(spec)})")
        if not (math.isfinite(phi) and 0.0 <= phi <= 1.0):
            raise KernelReferenceError(f"phi must be in [0, 1], got {phi!r}")
        if any(c for _name, c in spec.initial_mobile_counts):
            raise KernelReferenceError("initial mobile molecules are no state of the kernel model")
        unit = table["moment_unit"]
        self.ppu, self.nbb = int(unit["periods_per_declared_unit"]), int(unit["backbone_atoms_per_period"])
        self.unit_mass = float(unit["molar_mass_g_mol"])
        self.phi = float(phi)
        self.families = [f for f in K.FAMILIES if f in table["families"]]
        self.fam = {f: table["families"][f] for f in self.families}
        self.period_composition = {el: Fraction(v) for el, v in unit["composition"].items()}
        self.composition_offset = {f: {el: Fraction(v) for el, v in self.fam[f]["composition_offset"].items()}
                                   for f in self.families}
        self.species_of = {(f, int(key.rsplit("@", 1)[1])): name for f in self.families
                           for key, name in self.fam[f]["declared_species"].items()}
        init = declared_initialization(spec)
        initial = Counter()
        for dp, count in init.dp_counts:
            initial[dp * self.ppu] += count
        self.n_max = int(n_max or 4 * max(initial))
        if self.n_max < max(initial):
            raise KernelReferenceError(f"n_max {self.n_max} is below the longest initial chain ({max(initial)})")
        self.size = len(self.families) * self.n_max
        self.initial = numpy.zeros(self.size)
        for n, count in initial.items():
            self.initial[self.index("chain", n)] = count
        self.scale, self.volume = float(init.scale_factor), float(init.volume_m3)
        self.nav = constants.Na * self.volume
        self.u = self.scale / constants.Na              # mol per representative molecule
        self.a = spec.repeat_unit.composition.molar_mass_g_mol()
        self.e = spec.end_composition.molar_mass_g_mol()
        self.b_unit = self.nbb * self.ppu
        self._compile(table, kernels)

    # -- states ----------------------------------------------------------------
    def index(self, family, n) -> int:
        return self.families.index(family) * self.n_max + n - 1

    def state(self, i) -> Tuple[str, int]:
        return self.families[i // self.n_max], i % self.n_max + 1

    def label_of(self, family, n) -> Optional[str]:
        f = self.fam.get(family)
        if f is None or n < f["min_n"]:
            return None
        return f["labels"].get(str(n), f["label_far"])

    def molar_mass(self, family, n) -> float:
        return self.unit_mass * n + self.fam[family]["molar_mass_offset_g_mol"]

    def composition(self, family, n) -> Dict[str, Fraction]:
        """Atoms per element of state (family, n): n periods plus the family's offset."""
        offset = self.composition_offset[family]
        return {el: v for el in sorted(set(self.period_composition) | set(offset))
                if (v := n * self.period_composition.get(el, 0) + offset.get(el, 0))}

    def _target_composition(self, target) -> Counter:
        family, n = (self.families[target[1]], target[2]) if isinstance(target, tuple) else self.state(target)
        return Counter(self.composition(family, n))

    def _unbalanced(self, kernel, template, delta: Counter):
        delta = {el: v for el, v in delta.items() if v}
        if delta and not any(u["template"] == template for u in self.unbalanced):
            self.unbalanced.append({"kernel": kernel, "template": template,
                                    "delta": {el: str(v) for el, v in sorted(delta.items())}})

    def polymer(self, family, n) -> bool:
        label = self.label_of(family, n)
        return (label is not None and f"{label}@{n}" not in self.fam[family]["declared_species"]
                and (n * self.nbb) // self.b_unit >= 1)

    def _target(self, label, n):
        family = K.family_of(label)
        if family not in self.fam or self.label_of(family, n) != label:
            raise KernelReferenceError(f"a product {label} of {n} periods is no state of the table")
        return self.index(family, n) if n <= self.n_max else ("over", self.families.index(family), n)

    # -- laws ------------------------------------------------------------------
    def _outcomes(self, law, family, n) -> List[Tuple[float, Tuple]]:
        """[(multiplicity, products)] of one participant law on state (family, n):
        products are (label, n') tuples or a join tuple."""
        if self.label_of(family, n) != law["class"]:
            return []
        if law["anchor"] in ("chain", "end"):
            if n <= int(law["exact_to_n"]):
                return [(float(m), tuple(tuple(p) for p in prods))
                        for m, prods in law["exact"].get(str(n), [])]
            got = K.predict(law, n, self.nbb, lambda f, m: self.label_of(f, m) or f"{f}:none")
            return [(float(m), prods) for prods, m in sorted(got.items())]
        L = n * self.nbb
        xs = range(1, L - 1) if family == "radical_interior" else range(1, L - 2)
        out = Counter()
        for m, prods in law["far"]["outcomes"]:
            placed = {}
            for x in xs:
                coded = []
                for p in prods:
                    if p[0] == "join":
                        coded.append(tuple(p))
                        continue
                    f, alpha, beta, gamma = p
                    lp = alpha * L + beta * x + gamma
                    label = self.label_of(f, lp // self.nbb) if lp >= 1 and lp % self.nbb == 0 else None
                    if label is None:
                        break
                    coded.append((label, lp // self.nbb))
                else:
                    placed[x] = tuple(sorted(coded))
            for x, coded in placed.items():
                out[coded] += float(m) / len(placed)
        return [(m, prods) for prods, m in sorted(out.items())]

    def _compile(self, table, keep):
        self.uni: Dict[int, List[Tuple[float, Tuple]]] = defaultdict(list)   # state -> [(rate, targets)]
        self.bis: List[_Bi] = []
        self.pair_overrides: List[dict] = []
        self.not_modelled, self.inert, self.excluded, self.uni_terms = [], [], [], []
        self.unbalanced: List[dict] = []
        uni_overrides = []
        for kid, kernel in sorted(table["kernels"].items()):
            if keep is not None and kid not in keep:
                self.excluded.append(kid)
                continue
            for term in kernel["templates"]:
                name = term["template"]
                if term.get("law_inert"):
                    self.inert.append({"kernel": kid, "template": name, "reason": term["law_inert"]})
                    continue
                if term.get("law") is None:
                    self.not_modelled.append({"kernel": kid, "template": name,
                                              "reason": term.get("law_refused") or "no product law"})
                    continue
                phi = self.phi if term["phi_thinned"] else 1.0
                if kernel["molecularity"] == 1:
                    w = self._weight(table, name)
                    self._add_uni(term["law"], w * phi, kid, name)
                    self.uni_terms.append({"kernel": kid, "template": name})
                    uni_overrides += [(term["law"], w * phi, o, phi) for o in term.get("overrides", [])]
                else:
                    w = self._weight(table, name)
                    bi = self._add_bi(kid, name, term["law"], term["w_over_divisor"] * phi)
                    for o in term.get("overrides", []):
                        self.pair_overrides.append(self._pair_override(bi, o, w, term["w_over_divisor"] / w, phi))
        for law, w, o, phi in uni_overrides:
            self._apply_uni_override(law, w, o, phi)
        if keep is not None and set(keep) - set(table["kernels"]):
            raise KernelReferenceError(f"unknown kernels {sorted(set(keep) - set(table['kernels']))}")

    @staticmethod
    def _weight(table, name) -> float:
        row = next(r for r in table["templates"] if r["template"] == name)
        return float(row["weight"])

    def _add_uni(self, law, w, kernel=None, template=None):
        family = law["family"]
        for n in range(1, self.n_max + 1):
            for m, prods in self._outcomes(law, family, n):
                targets = tuple(self._target(label, n2) for label, n2 in prods)
                self.uni[self.index(family, n)].append((w * m, targets))
                delta = Counter()
                for tg in targets:
                    delta.update(self._target_composition(tg))
                delta.subtract(self.composition(family, n))
                self._unbalanced(kernel, template, delta)

    def _add_bi(self, kid, name, law, k) -> _Bi:
        laws = [law["S0"], law["S1"]]
        sites, outcomes, joins = [], [], set()
        for slot_law in laws:
            family = slot_law["family"]
            s = numpy.zeros(self.size)
            per = {}
            for n in range(1, self.n_max + 1):
                rows = self._outcomes(slot_law, family, n)
                if not rows:
                    continue
                i = self.index(family, n)
                s[i] = math.fsum(m for m, _p in rows)
                kept = []
                for m, prods in rows:
                    (p,) = prods
                    if p[0] == "join":
                        joins.add((p[1], int(p[2])))
                    else:
                        kept.append((m, self._target(*p)))
                per[i] = kept
            sites.append(s)
            outcomes.append(per)
        join = None
        if joins:
            if len(joins) != 1 or any(len(v) for per in outcomes for v in per.values()):
                raise KernelReferenceError(f"{name}: a join law that is not one product for every site pair")
            family, gamma = next(iter(joins))
            if gamma % self.nbb:
                raise KernelReferenceError(f"{name}: a joined product off the period grid")
            join = (self.families.index(family), gamma // self.nbb)
            # affine in n: n0 + n1 + g periods of the joined family against the two reactants
            delta = Counter(self.composition(family, gamma // self.nbb))
            for slot_law in laws:
                delta.subtract(self.composition(slot_law["family"], 0))
            self._unbalanced(kid, name, delta)
        else:
            changes = []
            for per in outcomes:
                seen = set()
                for i, kept in per.items():
                    for _m, tg in kept:
                        d = self._target_composition(tg)
                        d.subtract(self._target_composition(i))
                        seen.add(tuple(sorted((el, v) for el, v in d.items() if v)))
                changes.append(seen)
            if any(len(c) > 1 for c in changes):
                self._unbalanced(kid, name, Counter({"(participant change depends on the site)": 1}))
            else:
                total = Counter()
                for c in changes:
                    for el, v in next(iter(c), ()):
                        total[el] += v
                self._unbalanced(kid, name, total)
        bi = _Bi(kid, name, k, sites, outcomes, join)
        self.bis.append(bi)
        return bi

    def _state_of(self, reactant) -> List[Tuple[int, Optional[list]]]:
        label, n, positions = reactant
        if n > self.n_max:
            return None
        return self.index(K.family_of(label), n), positions

    def _apply_uni_override(self, law, w, o, phi):
        """Replace w by the override's weight on its matched transitions."""
        ((label, n, positions),) = o["reactants"]
        if n > self.n_max:
            return
        family = K.family_of(label)
        i = self.index(family, n)
        products = tuple(sorted((p[0], p[1]) for p in o["products"]))
        targets = tuple(self._target(*p) for p in products)
        w_o = o["weight"] * phi
        if law["anchor"] == "interior":
            L = n * self.nbb
            mirror = (lambda x: L - 1 - x) if family == "radical_interior" else (lambda x: L - 2 - x)
            xs = sorted({x for x in positions} | {mirror(x) for x in positions})
            share = 0.0
            for m, prods in law["far"]["outcomes"]:
                valid = {}
                for x in (range(1, L - 1) if family == "radical_interior" else range(1, L - 2)):
                    coded = []
                    for f, alpha, beta, gamma in prods:
                        lp = alpha * L + beta * x + gamma
                        lab = self.label_of(f, lp // self.nbb) if lp >= 1 and lp % self.nbb == 0 else None
                        if lab is None:
                            break
                        coded.append((lab, lp // self.nbb))
                    else:
                        valid[x] = tuple(sorted(coded))
                hits = [x for x in xs if valid.get(x) == products]
                if hits and o["multiplicity"] > m:
                    raise KernelReferenceError(f"override {o['name']}: more matched transitions than the law has")
                share += len(hits) * o["multiplicity"] / max(len(valid), 1)
            if not share:
                raise KernelReferenceError(f"override {o['name']}: its reaction is not one of the kernel model's")
            delta = (w_o - w) * share
        else:
            base = [m for m, prods in self._outcomes(law, family, n) if prods == products]
            if not base or o["multiplicity"] > base[0]:
                raise KernelReferenceError(f"override {o['name']}: its reaction is not one of the kernel model's")
            delta = (w_o - w) * o["multiplicity"]
        rows = self.uni[i]
        for k, (rate, tg) in enumerate(rows):
            if tg == targets:
                rows[k] = (rate + delta, tg)
                if rows[k][0] < -1e-12 * abs(rate):
                    raise KernelReferenceError(f"override {o['name']}: a negative rate")
                return
        raise KernelReferenceError(f"override {o['name']}: no channel of the kernel model makes its products")

    def _pair_override(self, bi: _Bi, o, w, inv_divisor, phi) -> dict:
        """A correction channel: (w_o - w)/divisor * multiplicity per
        unordered molecule pair of the two reactant states."""
        states = [self._state_of(r) for r in o["reactants"]]
        if any(s is None for s in states):
            return {"skip": True}
        if any(pos is not None for _i, pos in states):
            raise KernelReferenceError(f"override {o['name']}: a bimolecular override on an interior pool "
                                       "is not representable without positions")
        (a, _), (b, _) = states
        products = tuple(sorted((p[0], p[1]) for p in o["products"]))
        targets = tuple(self._target(*p) for p in products)
        # the kernel model's own rate of these transitions (per unordered pair): both slot assignments
        base = 0.0
        for i, j in ((a, b), (b, a)) if a != b else ((a, a),):
            if bi.join is not None:
                fj, g = bi.join
                n_out = self.state(i)[1] + self.state(j)[1] + g
                label = self.label_of(self.families[fj], n_out)
                if label is not None and (self._target(label, n_out),) == targets:
                    base += bi.sites[0][i] * bi.sites[1][j]
            else:
                for m0, t0 in bi.outcomes[0].get(i, []):
                    for m1, t1 in bi.outcomes[1].get(j, []):
                        if tuple(sorted((t0, t1), key=repr)) == tuple(sorted(targets, key=repr)):
                            base += m0 * m1
        if a == b:
            base *= 2.0                     # both slot assignments of an unordered pair of one state
        base *= bi.k
        delta = (o["weight"] - w) * inv_divisor * phi * o["multiplicity"]
        if base + delta < -1e-12 * base or base <= 0.0:
            raise KernelReferenceError(f"override {o['name']}: its reaction is not one of the kernel model's")
        return {"term": bi, "a": a, "b": b, "delta": delta, "base": base, "targets": targets,
                "name": o["name"]}

    # -- observables -----------------------------------------------------------
    def observables(self, counts: numpy.ndarray, over: Optional[numpy.ndarray] = None) -> dict:
        """whole_* (mol) of state counts of the representative system, with
        the overflow accumulators ((family, [N, sum n, sum n^2]))."""
        mu = [0.0, 0.0, 0.0]
        a, e = self.a, self.e
        for i in numpy.flatnonzero(counts):
            family, n = self.state(i)
            if not self.polymer(family, n):
                continue
            nm = (self.molar_mass(family, n) - e) / a
            c = counts[i]
            mu[0] += c
            mu[1] += c * nm
            mu[2] += c * nm * nm
        if over is not None:
            for f, family in enumerate(self.families):
                k0, k1, k2 = over[3 * f: 3 * f + 3]
                off = (self.fam[family]["molar_mass_offset_g_mol"] - e) / a
                s = self.unit_mass / a          # mass DP = s*n + off on this family
                mu[0] += k0
                mu[1] += s * k1 + off * k0
                mu[2] += s * s * k2 + 2 * s * off * k1 + off * off * k0
        mu = [self.u * x for x in mu]
        mn, mw = whole_averages(mu, a, e)
        return {"whole_mu0_mol": mu[0], "whole_mu1_mol": mu[1], "whole_mu2_mol": mu[2],
                "whole_Mn_g_mol": mn, "whole_Mw_g_mol": mw}

    def row_observation(self, counts: numpy.ndarray, over: Optional[numpy.ndarray] = None) -> dict:
        """The comparison-row inputs of state counts (representative system):
        atoms per element (exact integers when the counts are; else floats),
        mass (g, representative), declared species by name, polymer-derived
        states by length in periods (on the grid), overflow molecules, and
        the ``whole_*`` observables."""
        atoms = Counter()
        mass = []
        species = Counter()
        species_mass = []
        lengths = Counter()
        exact = True
        for i in numpy.flatnonzero(counts):
            family, n = self.state(i)
            c = counts[i]
            exact = exact and float(c).is_integer()
            k = int(c) if exact else Fraction(float(c))
            for el, v in self.composition(family, n).items():
                atoms[el] += k * v
            m = c * self.molar_mass(family, n)
            mass.append(m)
            name = self.species_of.get((family, n))
            if name is not None:
                species[name] += c
                species_mass.append(m)
            elif self.polymer(family, n):
                lengths[n] += c
        overflow = 0.0
        if over is not None:
            for f, family in enumerate(self.families):
                k0, k1, _k2 = over[3 * f: 3 * f + 3]
                if not k0:
                    continue
                exact = exact and float(k0).is_integer() and float(k1).is_integer()
                n0, n1 = (int(k0), int(k1)) if exact else (Fraction(float(k0)), Fraction(float(k1)))
                for el in set(self.period_composition) | set(self.composition_offset[family]):
                    atoms[el] += n1 * self.period_composition.get(el, 0) + n0 * self.composition_offset[family].get(el, 0)
                mass.append(self.unit_mass * k1 + self.fam[family]["molar_mass_offset_g_mol"] * k0)
                overflow += k0
        elements = {el: (int(v) if exact and Fraction(v).denominator == 1 else float(v))
                    for el, v in sorted(atoms.items())}
        return {"elements": elements, "exact_counts": exact, "mass_g": math.fsum(mass),
                "species": {k: float(v) for k, v in sorted(species.items())},
                "species_mass_g": math.fsum(species_mass),
                "lengths": {str(n): float(v) for n, v in sorted(lengths.items())},
                "overflow_molecules": float(overflow), "whole": self.observables(counts, over)}

    def mass_periods(self, counts, over=None) -> float:
        """Total backbone periods held (conserved by every kernel)."""
        n = numpy.tile(numpy.arange(1, self.n_max + 1, dtype=float), len(self.families))
        total = float(numpy.dot(n, counts))
        if over is not None:
            total += float(sum(over[3 * f + 1] for f in range(len(self.families))))
        return total


# --------------------------------------------------------------------------
# V0: the population balance
# --------------------------------------------------------------------------
class _Balance:
    """dN/dt and its Jacobian; y = [N (states), overflow (3 per family)]."""

    def __init__(self, model: KernelModel):
        self.m = model
        size, nf = model.size, len(model.families)
        self.size, self.nf = size, nf
        rows, cols, vals = [], [], []
        orow, ocol, oval = [], [], []

        def put(target, i, rate):
            if isinstance(target, tuple):
                _o, f, n = target
                for k, v in enumerate((1.0, n, n * n)):
                    orow.append(3 * f + k)
                    ocol.append(i)
                    oval.append(rate * v)
            else:
                rows.append(target)
                cols.append(i)
                vals.append(rate)

        for i, channels in model.uni.items():
            for rate, targets in channels:
                rows.append(i)
                cols.append(i)
                vals.append(-rate)
                for t in targets:
                    put(t, i, rate)
        self.A = sparse.csr_matrix((vals, (rows, cols)), shape=(size, size))
        self.Ao = sparse.csr_matrix((oval, (orow, ocol)), shape=(3 * nf, size))
        self.terms = []
        for bi in model.bis:
            gs, gos = [], []
            for slot in (0, 1):
                rows, cols, vals = [], [], []
                orow, ocol, oval = [], [], []
                for i, s in enumerate(bi.sites[slot]):
                    if s:
                        rows.append(i)
                        cols.append(i)
                        vals.append(-s)
                for i, rows_ in bi.outcomes[slot].items():
                    for mult, target in rows_:
                        put(target, i, mult)
                gs.append(sparse.csr_matrix((vals, (rows, cols)), shape=(size, size)))
                gos.append(sparse.csr_matrix((oval, (orow, ocol)), shape=(3 * nf, size)))
            self.terms.append((bi, gs, gos))

    def _join(self, bi, N):
        """Gain vector and overflow of a joined product, and per-slot length profiles."""
        m = self.m
        fj, g = bi.join
        u = [(bi.sites[s] * N).reshape(self.nf, m.n_max).sum(axis=0) for s in (0, 1)]
        conv = numpy.convolve(u[0], u[1])          # index k: n0 + n1 = k + 2
        gain = numpy.zeros(self.size)
        over = numpy.zeros(3 * self.nf)
        for k, v in enumerate(conv):
            if not v:
                continue
            n = k + 2 + g
            if n < 1:
                raise KernelReferenceError(f"{bi.template}: a joined product shorter than one period")
            if n <= m.n_max:
                gain[fj * m.n_max + n - 1] += v
            else:
                over[3 * fj: 3 * fj + 3] += v * numpy.array([1.0, n, n * n])
        return gain, over, u

    def rhs(self, _t, y):
        m = self.m
        N = y[:self.size]
        dN = self.A @ N
        do = self.Ao @ N
        for bi, gs, gos in self.terms:
            k = bi.k / m.nav
            S = [float(bi.sites[0] @ N), float(bi.sites[1] @ N)]
            dN += k * (S[1] * (gs[0] @ N) + S[0] * (gs[1] @ N))
            do += k * (S[1] * (gos[0] @ N) + S[0] * (gos[1] @ N))
            if bi.join is not None:
                gain, over, _u = self._join(bi, N)
                dN += k * gain
                do += k * over
        for po in m.pair_overrides:
            if po.get("skip"):
                continue
            rate = self._po_rate(po, N)
            dN[po["a"]] -= rate
            dN[po["b"]] -= rate
            for t in po["targets"]:
                if isinstance(t, tuple):
                    do[3 * t[1]: 3 * t[1] + 3] += rate * numpy.array([1.0, t[2], t[2] ** 2])
                else:
                    dN[t] += rate
        return numpy.concatenate([dN, do])

    def _po_rate(self, po, N):
        pairs = N[po["a"]] * N[po["b"]] if po["a"] != po["b"] else 0.5 * N[po["a"]] ** 2
        return po["delta"] * pairs / self.m.nav

    def jac(self, _t, y):
        m = self.m
        N = y[:self.size]
        J = numpy.zeros((self.size + 3 * self.nf, self.size + 3 * self.nf))
        J[:self.size, :self.size] = self.A.toarray()
        J[self.size:, :self.size] = self.Ao.toarray()
        for bi, gs, gos in self.terms:
            k = bi.k / m.nav
            S = [float(bi.sites[0] @ N), float(bi.sites[1] @ N)]
            for top, mats in ((slice(0, self.size), gs), (slice(self.size, None), gos)):
                J[top, :self.size] += k * (S[1] * mats[0].toarray() + S[0] * mats[1].toarray())
                J[top, :self.size] += k * (numpy.outer(mats[0] @ N, bi.sites[1])
                                           + numpy.outer(mats[1] @ N, bi.sites[0]))
            if bi.join is not None:
                fj, g = bi.join
                _gain, _over, u = self._join(bi, N)
                for slot in (0, 1):
                    other = u[1 - slot]
                    for i in numpy.flatnonzero(bi.sites[slot]):
                        n0 = i % m.n_max + 1
                        s = bi.sites[slot][i]
                        for n1 in numpy.flatnonzero(other) + 1:
                            n = n0 + n1 + g
                            v = k * s * other[n1 - 1]
                            if n <= m.n_max:
                                J[fj * m.n_max + n - 1, i] += v
                            else:
                                J[self.size + 3 * fj: self.size + 3 * fj + 3, i] += v * numpy.array([1.0, n, n * n])
        for po in m.pair_overrides:
            if po.get("skip"):
                continue
            a, b = po["a"], po["b"]
            d = po["delta"] / m.nav
            grads = {a: d * N[b], b: d * N[a]} if a != b else {a: d * N[a]}
            for col, gval in grads.items():
                J[a, col] -= gval
                J[b, col] -= gval
                for t in po["targets"]:
                    if isinstance(t, tuple):
                        J[self.size + 3 * t[1]: self.size + 3 * t[1] + 3, col] += gval * numpy.array([1.0, t[2], t[2] ** 2])
                    else:
                        J[t, col] += gval
        return J


def population_balance(model: KernelModel, times: Sequence[float], *, rtol: float = 1e-10,
                       atol: Optional[float] = None) -> List[dict]:
    """V0: the state counts at ``times`` (s) and their observables."""
    bal = _Balance(model)
    y0 = numpy.concatenate([model.initial, numpy.zeros(3 * len(model.families))])
    times = [float(t) for t in times]
    atol = atol if atol is not None else 1e-12 * float(model.initial.sum())
    if times[-1] > 0.0:
        sol = solve_ivp(bal.rhs, (0.0, times[-1]), y0, method="BDF", jac=bal.jac, t_eval=times,
                        rtol=rtol, atol=atol)
        if not sol.success:
            raise KernelReferenceError(f"the population balance failed: {sol.message}")
        ys = sol.y.T
    else:
        ys = [y0 for _t in times]
    out = []
    total = model.mass_periods(model.initial)
    for t, y in zip(times, ys):
        N, over = numpy.array(y[:model.size]), numpy.array(y[model.size:])
        row = {"time_s": t, **model.observables(N, over),
               "overflow_period_fraction": float(sum(over[3 * f + 1] for f in range(len(model.families))) / total),
               "counts": N, "overflow": over, "row": model.row_observation(N, over)}
        out.append(row)
    return out


# --------------------------------------------------------------------------
# V1: the stochastic simulation
# --------------------------------------------------------------------------
def law_domain_fractions(model: KernelModel, counts, overflow, validity: Sequence[dict]) -> List[dict]:
    """The law-domain number diagnostic of a MoM kernel set on one
    population: per declared term (the MoM producer's record ``validity``)
    and participant, the fraction of the participant family's chains (mu0;
    the overflow, longer than the grid, included) shorter than the length
    from which the term's moment law is exact (``participants``:
    ``exact_from_periods`` per slot; an older record's single
    ``valid_from_periods`` for every family).  None for a family the model
    lacks or that holds no chain; withheld (None, ``withheld``) when that
    length lies past the grid, where the overflow keeps no lengths."""
    out = []
    for v in validity:
        if v.get("participants"):
            slots = [(p["family"], int(p["exact_from_periods"]), {"slot": p["slot"]}) for p in v["participants"]]
        else:
            slots = [(family, int(v["valid_from_periods"]), {}) for family in v["families"]]
        for family, length, extra in slots:
            entry = {"kernel": v["kernel"], "template": v["template"], "family": family, **extra,
                     "valid_from_periods": int(v["valid_from_periods"])}
            if extra:
                entry["exact_from_periods"] = length
            total = below = None
            if family in model.families:
                j = model.families.index(family)
                block = numpy.asarray(counts[j * model.n_max:(j + 1) * model.n_max], dtype=float)
                total = float(block.sum() + overflow[3 * j])
                below = float(block[:max(0, min(length - 1, model.n_max))].sum())
                if length - 1 > model.n_max and overflow[3 * j] > 0.0:
                    entry["withheld"] = (f"the exact length {length} lies past the grid ({model.n_max} periods), "
                                         "where the overflow keeps no lengths")
                    total = None
            entry["mu0_fraction_below"] = below / total if total else None
            out.append(entry)
    return out


ATTRIBUTION_MOMENTS = (0, 1, 2)


def bind_mom_kernels(model: KernelModel, table: dict, record: dict) -> None:
    """Refuse, by name, a MoM producer's kernel record (``producer.json``
    facts.kernels.record) that is not of this table, temperature, geometry,
    phi and chemistry, names a template the model does not carry, or lacks
    the term laws the attribution reads."""
    if not isinstance(record, dict) or not isinstance(record.get("laws"), list) or "validity" not in record:
        raise KernelReferenceError("no MoM law record (facts.kernels.record.laws and validity of a MoM producer "
                                   "manifest): the attribution cannot be bound")
    if float(record.get("temperature_K", float("nan"))) != float(table["temperature_K"]):
        raise KernelReferenceError(f"the MoM record is at temperature {record.get('temperature_K')!r} K, the table "
                                   f"at {table['temperature_K']!r} K")
    mine, theirs = table["moment_unit"], record.get("moment_unit") or {}
    for key in ("periods_per_declared_unit", "backbone_atoms_per_period", "molar_mass_g_mol"):
        if theirs.get(key) != mine[key]:
            raise KernelReferenceError(f"the MoM record's period geometry {key} = {theirs.get(key)!r} is not the "
                                       f"table's {mine[key]!r}")
    chemistry = (table.get("catalog") or {}).get("chemistry_sha256")
    if record.get("chemistry_sha256") != chemistry:
        raise KernelReferenceError(f"the MoM record's chemistry {record.get('chemistry_sha256')!r} is not the "
                                   f"table's {chemistry!r}")
    if float(record.get("phi", float("nan"))) != model.phi:
        raise KernelReferenceError(f"the MoM record's phi {record.get('phi')!r} is not the reference's {model.phi!r}")
    carried = {t["template"]: kid for kid, k in table["kernels"].items() if kid not in model.excluded
               for t in k["templates"]}
    names = [n for law in record["laws"] for n in law["templates"]]
    missing = sorted({n for n in names if n not in carried})
    if missing:
        raise KernelReferenceError(f"the MoM record's templates {missing} are not terms the reference carries")
    if K.template_rows_sha256(table, names) != record.get("table_templates_sha256"):
        raise KernelReferenceError("the MoM record's template rows (weights, sites, laws) are not this table's "
                                   f"(table_templates_sha256 {record.get('table_templates_sha256')!r})")
    terms = {(v["kernel"], v["template"]) for v in record["validity"]}
    unlisted = sorted(law["template"] for law in record["laws"] if (law["kernel"], law["template"]) not in terms)
    if unlisted:
        raise KernelReferenceError(f"the MoM record's terms {unlisted} have no validity entry (a partial record)")


def read_mom_kernels(path) -> dict:
    """The MoM kernel record of a producer manifest (``producer.json``)."""
    doc = json.loads(Path(path).read_text())
    record = (((doc.get("facts") or {}).get("kernels") or {}).get("record") if isinstance(doc, dict) else None)
    if not isinstance(record, dict):
        raise KernelReferenceError(f"{path}: not a MoM producer manifest with facts.kernels.record")
    return record


def _site_series(model: KernelModel, law: dict, family: str, product: Optional[str]) -> Tuple[numpy.ndarray, Optional[tuple]]:
    """A participant law's sites on (family, n) as the balance reads them
    (``_outcomes``), counting only outcomes into ``product`` when given: the
    grid's values and, past the grid, (gamma, delta) with sites = gamma +
    delta n (None when they are not affine there)."""
    key = (id(law), family, product)
    cache = model.__dict__.setdefault("_site_series", {})
    if key in cache:
        return cache[key]

    def at(n):
        total = 0.0
        for m, prods in model._outcomes(law, family, n):
            if product is None or (len(prods) == 1 and prods[0][0] != "join"
                                   and K.family_of(prods[0][0]) == product):
                total += m
        return total
    grid = numpy.array([at(n) for n in range(1, model.n_max + 1)])
    a, b, c = (at(model.n_max + i) for i in (1, 2, 3))
    delta = b - a
    beyond = (a - delta * (model.n_max + 1), delta) if abs(c - b - delta) <= 1e-12 * max(abs(a), abs(b), abs(c), 1.0) \
        else None
    cache[key] = (grid, beyond)
    return grid, beyond


def _weighted(model: KernelModel, family: str, grid: numpy.ndarray, beyond, k: int, counts, overflow,
              o3: Optional[float] = None):
    """(sum over the family's chains of f(n) n^k / ppu^k, why not): the grid
    exactly, the overflow from its mu0..mu2 where f is affine past the grid
    (its sum of n^3, which the overflow does not keep, given as ``o3``)."""
    j = model.families.index(family)
    n = numpy.arange(1, model.n_max + 1, dtype=float)
    block = numpy.asarray(counts[j * model.n_max:(j + 1) * model.n_max], dtype=float)
    total = float(numpy.dot(grid * n ** k, block)) / model.ppu ** k
    o = [float(x) for x in overflow[3 * j:3 * j + 3]]
    if o[0] == 0.0:
        return total, None
    if beyond is None:
        return None, f"{o[0]!r} chains of {family} past the grid, where its sites are not affine"
    gamma, delta = beyond
    if delta != 0.0 and k + 1 > 2 and o3 is None:
        return None, f"{o[0]!r} chains of {family} past the grid: n^{k + 1} needs the overflow's third moment"
    over = gamma * o[k] + (delta * (o[k + 1] if k + 1 <= 2 else o3) if delta != 0.0 else 0.0)
    return total + over / model.ppu ** k, None


def law_attribution(model: KernelModel, table: dict, record: dict, counts, overflow) -> List[dict]:
    """The law-domain attribution of a MoM kernel set on one population of
    the balance: per MoM term (``record["laws"]``: the moment law MoM
    integrates, per chain or per pair in the balance's units) and
    participant, the signed difference between the law's propensity and the
    exact one (the balance's own sites of the term's templates, into the
    term's product families), summed over the population and weighted by
    the participant's length^k (declared units), k = 0, 1, 2 -- the rate at
    which the term moves mu_k out of the participant's family -- in mol/s;
    with the exact value and the number fraction below the participant's
    exact-from length.  A moment the overflow cannot carry is withheld
    (None, ``withheld``)."""
    entries = {t["template"]: t for k in table["kernels"].values() for t in k["templates"]}
    rows = {r["template"]: r for r in table["templates"]}
    fractions = {(f["template"], f.get("slot")): f for f in law_domain_fractions(model, counts, overflow,
                                                                                    record["validity"])}
    # the overflow's sum of n^3 per family, bounded: at least o2^2 / o1 (Cauchy-Schwarz) and (n_max + 1) o2,
    # at most (all periods of the system) o2 -- no chain is longer than every period there is
    n = numpy.arange(1, model.n_max + 1, dtype=float)
    periods = sum(float(numpy.dot(n, counts[j * model.n_max:(j + 1) * model.n_max])) + float(overflow[3 * j + 1])
                  for j in range(len(model.families)))
    o3 = {}
    for j, family in enumerate(model.families):
        o0, o1, o2 = (float(x) for x in overflow[3 * j:3 * j + 3])
        o3[family] = (0.0, 0.0) if o0 == 0.0 else (max(o2 * o2 / o1, (model.n_max + 1) * o2), periods * o2)
    one = (numpy.ones(model.n_max), (1.0, 0.0))
    ramp = (numpy.arange(1, model.n_max + 1, dtype=float), (0.0, 1.0))
    out = []
    for law in record["laws"]:
        fams, products, c = law["families"], law.get("products"), law["coefficients"]
        phi = [model.phi if entries[t]["phi_thinned"] else 1.0 for t in law["templates"]]
        if law["kind"] == "uni":
            (family,) = fams
            product = products[0] if products else None
            exact = [(float(rows[t]["weight"]) * f, [_site_series(model, entries[t]["law"], family, product)])
                     for t, f in zip(law["templates"], phi)]
            law_parts = [(c[0], [one]), (c[1], [ramp])]
            scale = model.u
        else:
            exact = []
            for t, f in zip(law["templates"], phi):
                slot_laws = (entries[t]["law"]["S0"], entries[t]["law"]["S1"])
                exact.append((float(entries[t]["w_over_divisor"]) * f,
                              [_site_series(model, slot_laws[p], fams[p], products[p] if products else None)
                               for p in (0, 1)]))
            law_parts = [(c[0], [one, one]), (c[1], [ramp, one]), (c[2], [one, ramp]), (c[3], [ramp, ramp])]
            scale = model.u / model.nav
        slots = []
        for p, family in enumerate(fams):
            moments, withheld = {}, {}
            for k in ATTRIBUTION_MOMENTS:
                def total(parts, third):
                    acc = 0.0
                    for coef, series in parts:
                        if coef == 0.0:
                            continue
                        prod = coef
                        for q, (grid, beyond) in enumerate(series):
                            value, why = _weighted(model, fams[q], grid, beyond, k if q == p else 0, counts,
                                                   overflow, third)
                            if value is None:
                                return None, why
                            prod *= value
                        acc += prod
                    return acc * scale, None
                # linear in the weighted family's overflow sum of n^3: evaluated at both of its bounds
                (x_lo, why_x), (m_lo, why_m) = total(exact, o3[family][0]), total(law_parts, o3[family][0])
                (x_hi, _w), (m_hi, _v) = total(exact, o3[family][1]), total(law_parts, o3[family][1])
                if x_lo is None or m_lo is None:
                    withheld[str(k)] = why_x or why_m
                    moments[str(k)] = {"law_minus_exact_mol_s": None, "exact_mol_s": None}
                elif abs(x_hi - x_lo) <= 1e-12 * abs(x_lo) and abs((m_hi - x_hi) - (m_lo - x_lo)) <= 1e-12 * max(
                        abs(m_lo - x_lo), abs(x_lo)):
                    moments[str(k)] = {"law_minus_exact_mol_s": m_lo - x_lo, "exact_mol_s": x_lo}
                else:
                    moments[str(k)] = {"law_minus_exact_mol_s": None, "exact_mol_s": None,
                                       "law_minus_exact_bounds_mol_s": sorted([m_lo - x_lo, m_hi - x_hi]),
                                       "exact_bounds_mol_s": sorted([x_lo, x_hi])}
                    withheld[str(k)] = (f"chains of {family} past the grid: their sum of n^3 (the overflow keeps "
                                        f"n^0..n^2) only bounded, {list(o3[family])}; the value lies in the bounds")
            f = fractions.get((law["template"], p)) or {}
            slots.append({"slot": p, "family": family, "exact_from_periods": f.get("exact_from_periods"),
                          "mu0_fraction_below": f.get("mu0_fraction_below"), "moments": moments,
                          "withheld": withheld})
        out.append({"kernel": law["kernel"], "template": law["template"], "kind": law["kind"], "slots": slots})
    return out


def read_validity(path) -> List[dict]:
    """The validity lengths of a MoM producer manifest (``producer.json``:
    facts.kernels.record.validity) or of a bare JSON list of them."""
    doc = json.loads(Path(path).read_text())
    if isinstance(doc, dict):
        doc = (((doc.get("facts") or {}).get("kernels") or {}).get("record") or {}).get("validity")
    if not isinstance(doc, list) or not all(
            isinstance(v, dict) and {"kernel", "template", "families", "valid_from_periods"} <= set(v) for v in doc):
        raise KernelReferenceError(f"{path}: no MoM kernel validity lengths (a producer.json with "
                                   "facts.kernels.record.validity, or a list of them)")
    return doc


def _positive_overrides(model) -> list:
    return [po for po in model.pair_overrides if not po.get("skip") and po["delta"] > 0]


def propensities(model: KernelModel, N) -> numpy.ndarray:
    """The SSA channel propensities (1/s) of state counts ``N``: every
    unimolecular channel together, each bimolecular term
    ``k*(S0*S1 - X)/(N_A*V)``, each override channel of a higher weight."""
    out = [float(sum(math.fsum(r for r, _t in rows) * N[i] for i, rows in model.uni.items()))]
    for bi in model.bis:
        s0, s1 = bi.sites
        X = float((s0 * s1) @ N)
        out.append(max(bi.k * (float(s0 @ N) * float(s1 @ N) - X) / model.nav, 0.0))
    for po in _positive_overrides(model):
        pairs = N[po["a"]] * N[po["b"]] if po["a"] != po["b"] else 0.5 * N[po["a"]] * (N[po["a"]] - 1)
        out.append(po["delta"] * pairs / model.nav)
    return numpy.array(out)


def stochastic(model: KernelModel, times: Sequence[float], seed: int, *,
               max_events: Optional[int] = None) -> Tuple[List[dict], int]:
    """V1: one Gillespie run of the representative system; (rows at
    ``times``, events fired)."""
    rng = numpy.random.default_rng(seed)
    N = numpy.array(model.initial, dtype=float)
    over = numpy.zeros(3 * len(model.families))
    size = model.size
    K_uni = numpy.zeros(size)
    for i, rows in model.uni.items():
        K_uni[i] = math.fsum(r for r, _t in rows)
    times = [float(t) for t in times]
    out, t, events, nxt = [], 0.0, 0, 0
    pos = _positive_overrides(model)
    thin = defaultdict(list)
    for po in model.pair_overrides:
        if not po.get("skip") and po["delta"] < 0:
            thin[id(po["term"])].append(po)

    def record():
        out.append({"time_s": times[nxt], **model.observables(N, over), "row": model.row_observation(N, over)})

    def make(target):
        if isinstance(target, tuple):
            _o, f, n = target
            over[3 * f: 3 * f + 3] += (1.0, n, n * n)
        else:
            N[target] += 1

    def draw(weights):
        c = numpy.cumsum(weights)
        return int(numpy.searchsorted(c, rng.random() * c[-1], side="right"))

    while nxt < len(times):
        props = numpy.concatenate([[float(K_uni @ N)], propensities(model, N)[1:]])
        total = float(props.sum())
        dt = rng.exponential(1.0 / total) if total > 0 else math.inf
        while nxt < len(times) and t + dt > times[nxt]:
            record()
            nxt += 1
        if nxt >= len(times):
            break
        t += dt
        if max_events is not None and events >= max_events:
            raise KernelReferenceError(f"the stochastic run reached {max_events} events before the horizon")
        channel = draw(props)
        events += 1
        if channel == 0:
            i = draw(K_uni * N)
            rows = model.uni[i]
            rate, targets = rows[draw(numpy.array([r for r, _t in rows]))]
            N[i] -= 1
            for tg in targets:
                make(tg)
        elif channel <= len(model.bis):
            bi = model.bis[channel - 1]
            s0, s1 = bi.sites
            while True:
                i, j = draw(s0 * N), draw(s1 * N)
                if i != j or rng.random() * N[i] < N[i] - 1:
                    break
            if bi.join is None:
                r0 = bi.outcomes[0][i]
                r1 = bi.outcomes[1][j]
                t0 = r0[draw(numpy.array([m for m, _t in r0]))][1]
                t1 = r1[draw(numpy.array([m for m, _t in r1]))][1]
                made = (t0, t1)
            else:
                fj, g = bi.join
                n = model.state(i)[1] + model.state(j)[1] + g
                made = (fj * model.n_max + n - 1 if n <= model.n_max else ("over", fj, n),)
            for po in thin.get(id(bi), ()):
                if {i, j} == {po["a"], po["b"]} and sorted(made, key=repr) == sorted(po["targets"], key=repr):
                    if rng.random() * po["base"] >= po["base"] + po["delta"]:
                        made = None
                    break
            if made is None:
                events -= 1                         # thinned: the override's lower weight
                continue
            N[i] -= 1
            N[j] -= 1
            for tg in made:
                make(tg)
        else:
            po = pos[channel - 1 - len(model.bis)]
            N[po["a"]] -= 1
            N[po["b"]] -= 1
            for tg in po["targets"]:
                make(tg)
    return out, events


def ensemble(model: KernelModel, times, seeds: Sequence[int], **kw) -> dict:
    """Mean and standard error of every whole_* key over seeds."""
    runs = [stochastic(model, times, s, **kw) for s in seeds]
    rows = []
    for k, t in enumerate(times):
        row = {"time_s": float(t)}
        for key in WHOLE_ENSEMBLE_KEYS:
            values = numpy.array([r[0][k][key] for r in runs])
            se = float(values.std(ddof=1) / math.sqrt(len(values))) if len(values) > 1 else float("nan")
            row[key] = {"mean": float(values.mean()), "se": se}
        rows.append(row)
    return {"seeds": list(seeds), "events": [r[1] for r in runs], "rows": rows,
            "row_observations": [[x["row"] for x in r[0]] for r in runs]}


def comparison(pbe_rows, ssa) -> List[dict]:
    """Per time and whole_* key: the PBE value, the SSA mean and standard
    error, and z = (mean - pbe) / se (0 when both are exact and equal)."""
    out = []
    for p, s in zip(pbe_rows, ssa["rows"]):
        for key in WHOLE_ENSEMBLE_KEYS:
            mean, se = s[key]["mean"], s[key]["se"]
            diff = mean - p[key]
            scale = max(abs(p[key]), abs(mean), 1e-300)
            if se > 1e-12 * scale:
                z = diff / se
            else:
                z = 0.0 if abs(diff) <= 1e-9 * scale else math.inf
            out.append({"time_s": p["time_s"], "quantity": key, "pbe": p[key], "ssa_mean": mean,
                        "ssa_se": se, "z": z})
    return out


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------
def main(argv=None) -> int:
    from rmgpy.kmc import run as kr
    parser = argparse.ArgumentParser(prog="python -m rmgpy.kmc.kernel_reference",
                                     description="exact population balance and kernel SSA of a kernel-rates table")
    parser.add_argument("deck")
    parser.add_argument("table", help=f"a {K.SCHEMA} JSON (python -m rmgpy.kmc.kernel_rates)")
    parser.add_argument("--end", type=float, required=True)
    parser.add_argument("--snapshots", type=int, required=True)
    parser.add_argument("--seeds", type=int, default=8)
    parser.add_argument("--phi", type=float, default=1.0, help="cage-escape factor of radical-creating terms")
    parser.add_argument("--n-max", type=int, default=None, help="longest state of the PBE grid (periods)")
    parser.add_argument("--kernel", action="append", default=None, help="keep only these kernels (repeatable)")
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--out", default=None, help="default: <table stem>.reference.json beside it")
    parser.add_argument("--allow-unbound-table", action="store_true")
    parser.add_argument("--validity", default=None,
                        help="a MoM producer.json: per MoM term, the law-domain attribution (law minus exact "
                             "propensity, weighted by mu0..mu2) and the fraction of each participant family's "
                             "chains below its exact-from length (law_domain), bound to the producer's record")
    args = parser.parse_args(argv)
    try:
        deck = kr.read_deck(args.deck)
        text = Path(args.table).read_text()
        table = json.loads(text)
        model = KernelModel(table, deck.spec, phi=args.phi, n_max=args.n_max, kernels=args.kernel,
                            allow_unbound=args.allow_unbound_table)
        times = kr.observation_times(None, args.end, args.snapshots)
        record = read_mom_kernels(args.validity) if args.validity else None
        if record is not None:
            bind_mom_kernels(model, table, record)
        pbe = population_balance(model, times)
        ssa = ensemble(model, times, list(range(args.seeds)), max_events=args.max_events)
    except (KernelReferenceError, K.KernelRatesError, kr.RunError, OSError, ValueError) as error:
        print(f"refused: {error}", file=sys.stderr)
        return 2
    from rmgpy.kmc.melt_engine import SpeciesIdentity
    from rmgpy.kmc.melt_observables import declared_species_names
    observations = ssa.pop("row_observations")
    doc = {"schema": SCHEMA, "table": str(Path(args.table).resolve()),
           "table_sha256": hashlib.sha256(text.encode()).hexdigest(),
           "deck": {"path": str(Path(args.deck).resolve()), "sha256": deck.sha256,
                    "chemistry_sha256": kr.chemistry_sha256(deck.spec), "cutoff": deck.xs},
           "temperature_K": table["temperature_K"], "times_s": list(times), "phi": model.phi,
           "pass": "well_mixed_bound" if model.phi == 1.0 else f"phi={model.phi!r}", "n_max": model.n_max,
           "volume_m3": model.volume, "scale_factor": model.scale, "unit_mol": model.u,
           "geometry": {"periods_per_declared_unit": model.ppu, "backbone_atoms_per_period": model.nbb},
           "declared_species": declared_species_names(deck.spec, SpeciesIdentity(deck.spec)),
           "excluded_kernels": model.excluded, "inert_terms": model.inert, "not_modelled_terms": model.not_modelled,
           "unbalanced_terms": model.unbalanced,
           "row_observations": {"pbe": [r["row"] for r in pbe], "ssa": {"seeds": ssa["seeds"], "runs": observations}},
           "pbe": [{k: v for k, v in r.items() if k not in ("counts", "overflow", "row")} for r in pbe],
           "pbe_final_distribution": {f"{model.label_of(*model.state(i))}@{model.state(i)[1]}": float(v)
                                      for i, v in enumerate(pbe[-1]["counts"]) if v},
           "ssa": ssa, "comparison": comparison(pbe, ssa)}
    if record is not None:
        doc["law_domain"] = {
            "source": str(Path(args.validity).resolve()),
            "identity": {k: record.get(k) for k in ("table_templates_sha256", "chemistry_sha256", "temperature_K",
                                                    "moment_unit", "phi", "table_sha256")},
            "validity": record["validity"],
            "pbe": [law_domain_fractions(model, r["counts"], r["overflow"], record["validity"]) for r in pbe],
            "attribution": [law_attribution(model, table, record, r["counts"], r["overflow"]) for r in pbe]}
    out = Path(args.out) if args.out else Path(args.table).with_suffix(".reference.json")
    out.write_text(json.dumps(doc, indent=1, sort_keys=True, default=float) + "\n")
    worst = max((abs(r["z"]) for r in doc["comparison"]), default=0.0)
    print(f"written: {out} (max |z| {worst:.3g} over {len(doc['comparison'])} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
